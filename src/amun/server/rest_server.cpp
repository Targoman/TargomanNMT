/******************************************************************************
 * Targoman: A robust Machine Translation framework               *
 *                                                                            *
 * Copyright 2014-2018 by ITRC <http://itrc.ac.ir>                            *
 *                                                                            *
 * This file is part of Targoman.                                             *
 *                                                                            *
 * Targoman is free software: you can redistribute it and/or modify           *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Targoman is distributed in the hope that it will be useful,                *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with Targoman. If not, see <http://www.gnu.org/licenses/>.           *
 *                                                                            *
 ******************************************************************************/
/**
 * @author Behrooz Vedadian <vedadian@targoman.com>
 */
 
#include <stdexcept>
#include <iostream>
#include <memory>
#include <chrono>
#include <thread>
#include <atomic>
#include <vector>
#include <algorithm>
#include <string>
#include <string.h>
#include <regex>
#include <functional>
#include <evhttp.h>

#include <cstdlib>
#include <boost/timer/timer.hpp>

#include "common/god.h"
#include "common/logging.h"
#include "common/search.h"
#include "common/threadpool.h"
#include "common/printer.h"
#include "common/sentence.h"
#include "common/sentences.h"
#include "common/exception.h"
#include "common/translation_task.h"

#include <sys/socket.h>
#include <arpa/inet.h>

#include "gason.h"

using namespace amunmt;

static std::string SERVER_NAME = "";

void GetPrimaryIp(char* buffer, size_t buflen) 
{
    assert(buflen >= 16);

    int sock = socket(AF_INET, SOCK_DGRAM, 0);
    assert(sock != -1);

    const char* kGoogleDnsIp = "74.125.115.104";
    uint16_t kDnsPort = 80;
    struct sockaddr_in serv;
    memset(&serv, 0, sizeof(serv));
    serv.sin_family = AF_INET;
    serv.sin_addr.s_addr = inet_addr(kGoogleDnsIp);
    serv.sin_port = htons(kDnsPort);

    int err = connect(sock, (const sockaddr*) &serv, sizeof(serv));
    assert(err != -1);

    sockaddr_in name;
    socklen_t namelen = sizeof(name);
    err = getsockname(sock, (sockaddr*) &name, &namelen);
    assert(err != -1);

    const char* p = inet_ntop(AF_INET, &name.sin_addr, buffer, buflen);
    assert(p);

    close(sock);
}

std::string getPostBody(evhttp_request *_request)
{
    std::string result;
    evbuffer *buf = evhttp_request_get_input_buffer(_request);
    size_t len = evbuffer_get_length(buf);
    result.resize(len);
    evbuffer_copyout(buf, (void*)result.data(), len);
    return result;
}

typedef std::tuple<std::vector<std::string>, std::vector<std::string>> InputLines_t;

InputLines_t getTextLines(std::string _body)
{
    std::vector<std::string> srcs;
    std::vector<std::string> origs;
    Json::JsonValue value;
    Json::JsonAllocator allocator;
    char *s = (char *)_body.data();
    char *endptr;
    int status = Json::jsonParse(s, &endptr, &value, allocator);
    if (status != Json::JSON_OK)
        throw std::runtime_error("Input is not a valid JSON string.");
    if(value.getTag() != Json::JSON_ARRAY)
        throw std::runtime_error("Input string must be a JSON array.");
    for(auto e : value) {
        if(e->value.getTag() != Json::JSON_OBJECT)
            throw std::runtime_error("Each element of the input array must be an object containing `src` and `orig` fields.");
        char *src = nullptr, *orig = nullptr;
        for(auto i : e->value) {
            if(strcmp(i->key, "src") == 0) {
                if(i->value.getTag() != Json::JSON_STRING)
                    throw std::runtime_error("`src` field must contain string.");
                src = i->value.toString();
            }
            // if(strcmp(i->key, "orig") == 0) {
            //     if(i->value.getTag() != Json::JSON_STRING)
            //         throw std::runtime_error("`orig` field must contain string.");
            //     orig = i->value.toString();
            // }
        }
        srcs.push_back(src);
        if(orig == nullptr)
            origs.push_back(std::string());
        else
            origs.push_back(orig);
    }
    return InputLines_t(srcs, origs);
}

template<typename T>
void evbuffer_add_printf_array(evbuffer* _buff, char* _format, const T& _array) {
    evbuffer_add_printf(_buff, "[");
    if(_array.size() > 0) {
        evbuffer_add_printf_array(_buff, _format, _array[0]);
        for(size_t Index = 1; Index < _array.size(); ++Index) {
            evbuffer_add_printf(_buff, ",");
            evbuffer_add_printf_array(_buff, _format, _array[Index]);
        }
    }
    evbuffer_add_printf(_buff, "]");
}

template<>
void evbuffer_add_printf_array<std::string>(evbuffer* _buff, char* _format, const std::string& _array) {
    evbuffer_add_printf(_buff, _format, (char *)_array.data());
}

template<>
void evbuffer_add_printf_array<int>(evbuffer* _buff, char* _format, const int& _array) {
    evbuffer_add_printf(_buff, _format, _array);
}

void sendResults(
    evhttp_request *_request,
    std::vector<std::vector<std::string>> _source,
    std::vector<std::vector<std::vector<std::string>>> _phrases,
    std::vector<std::vector<std::vector<int>>> _alignments
) {
    auto *OutBuf = evhttp_request_get_output_buffer(_request);
    if (!OutBuf) {
        evhttp_send_error(_request, HTTP_INTERNAL, "{ \"code\": 500, \"msg\": \"Internal server error.\" }");
        return;        
    }
    
    bool First = true;
    evbuffer_add_printf(OutBuf, "{\"rslt\":[");
    for(size_t Index = 0; Index < _source.size(); ++Index) {
        if(!First)
            evbuffer_add_printf(OutBuf, ",");
        std::vector<std::string>& Tokens = _source[Index];
        std::vector<std::vector<std::string>>& Phrases = _phrases[Index];
        std::vector<std::vector<int>> Alignments = _alignments[Index];
        evbuffer_add_printf(OutBuf, "{");
        evbuffer_add_printf(OutBuf, "\"tokens\":");
        evbuffer_add_printf_array(OutBuf, "\"%s\"", Tokens);
        evbuffer_add_printf(OutBuf, ",\"phrases\":");
        evbuffer_add_printf_array(OutBuf, "\"%s\"", Phrases);
        evbuffer_add_printf(OutBuf, ",\"alignments\":");
        evbuffer_add_printf_array(OutBuf, "%d", Alignments);
        evbuffer_add_printf(OutBuf, "}");
        First = false;
    }
    evbuffer_add_printf(OutBuf, "],");
    evbuffer_add_printf(OutBuf, "\"serverName\":\"%s\"}", SERVER_NAME.c_str());
    evhttp_send_reply(_request, HTTP_OK, "", OutBuf);
}

typedef std::vector<std::vector<std::string>> NBestTranslations_t;
struct Alignment_t {
    size_t SourceSize;
    size_t TargetSize;
    std::vector<float> M;
};

void saveAlignmentsToList(const HypothesisPtr& hypothesis, size_t target_size, size_t source_size, Alignment_t& alignment) {
  std::vector<SoftAlignment> aligns;
  HypothesisPtr last = hypothesis->GetPrevHyp();
  while (last->GetPrevHyp().get() != nullptr) {
    SoftAlignment copy = *(last->GetAlignment(0));
    copy.resize(source_size - 1);
    aligns.insert(aligns.begin(), copy);
    last = last->GetPrevHyp();
  }
  alignment.SourceSize = source_size - 1;
  alignment.TargetSize = target_size - 1;
  alignment.M.resize((source_size - 1) * (target_size - 1));
  for(int i = 0; i < (int)target_size - 1; ++i) {
    const SoftAlignment& item = aligns[i];
    for(int j = 0; j < (int)item.size(); ++j) {
      const float value = item[j];
      alignment.M[i * (source_size - 1) + j] = value;
    }
  }
}

void saveHistoryToOutput(
  const History& history,
  const Sentence& sentence,
  God& __GOD__,
  std::vector<std::vector<std::string>>& _sentenceWords,
  std::vector<NBestTranslations_t>& _translations,
  std::vector<std::vector<Alignment_t>>& _alignments
) {
  God& god = __GOD__;
//   const Words& SentenceWordIds = sentence.GetWords(0);
//   const Vocab& SourceVocab = god.GetSourceVocab(0);
//   std::vector<std::string> ThisSentenceWords;
//   ThisSentenceWords.resize(SentenceWordIds.size() - 1);
//   for(size_t WordIndex = 0; WordIndex < ThisSentenceWords.size(); ++WordIndex)
//     ThisSentenceWords[WordIndex] = SourceVocab[SentenceWordIds[WordIndex]];
//   _sentenceWords.push_back(ThisSentenceWords);
  if (god.Get<bool>("n-best")) {
    std::vector<std::string> scorerNames = god.GetScorerNames();
    const NBestList &nbl = history.NBest(god.Get<size_t>("beam-size"));

    NBestTranslations_t sentence_outputs;
    std::vector<Alignment_t> sentence_alignment;

    for (size_t i = 0; i < nbl.size(); ++i) {
      const Result& result = nbl[i];
      const Words &words = result.first;
      const HypothesisPtr &hypo = result.second;

      std::vector<std::string> translation = god.Postprocess(god.GetTargetVocab()(words));
      sentence_outputs.push_back(translation);
      Alignment_t alignment;
      saveAlignmentsToList(hypo, translation.size() + 1, sentence.size(), alignment);
      sentence_alignment.push_back(alignment);
    }
    _translations.push_back(sentence_outputs);
    _alignments.push_back(sentence_alignment);
  } else {
    auto bestTranslation = history.Top();
    std::vector<std::string> bestSentenceWords = god.Postprocess(god.GetTargetVocab()(bestTranslation.first));

    NBestTranslations_t sentence_outputs;
    std::vector<Alignment_t> sentence_alignment;

    sentence_outputs.push_back(bestSentenceWords);
    Alignment_t alignment;
    saveAlignmentsToList(bestTranslation.second, bestSentenceWords.size() + 1, sentence.size(), alignment);
    sentence_alignment.push_back(alignment);

    _translations.push_back(sentence_outputs);
    _alignments.push_back(sentence_alignment);
  }
}

void translate(std::vector<std::string>& _lines, God& __GOD__, std::vector<std::vector<std::string>>& _sentenceWords, std::vector<NBestTranslations_t>& _translations, std::vector<std::vector<Alignment_t>>& _alignments)
{
//   LOG(info)->info("Retrieving batching settings.");

  size_t miniSize = __GOD__.Get<size_t>("mini-batch");
  size_t maxiSize = __GOD__.Get<size_t>("maxi-batch");
  int miniWords = __GOD__.Get<int>("mini-batch-words");

  std::vector<std::future< std::shared_ptr<Histories> >> results;
  SentencesPtr maxiBatch(new Sentences());
  std::vector<SentencePtr> sentences_by_lineno;

//   LOG(info)->info("Creating batches.");

  sentences_by_lineno.resize(_lines.size());
  for(int lineNum = 0; lineNum < _lines.size(); ++lineNum) {
    std::string& line = _lines[lineNum];

    std::vector<std::string> __words;
    Split(line, __words, " ");
    for(std::string& w : __words) {
        Word __index = __GOD__.GetSourceVocab(0)[w];
        if(__index == UNK_ID) {
            std::string __lower = w;
            std::transform(__lower.begin(), __lower.end(), __lower.begin(), ::tolower);
            __index = __GOD__.GetSourceVocab(0)[__lower];
            if(__index != UNK_ID)
                w = __lower;
        }
    }
    _sentenceWords.push_back(__words);
    SentencePtr __s(new Sentence(__GOD__, lineNum, line));
    sentences_by_lineno[lineNum] = __s;
    maxiBatch->push_back(__s);

    if (maxiBatch->size() >= maxiSize) {

      maxiBatch->SortByLength();
      while (maxiBatch->size()) {
        SentencesPtr miniBatch = maxiBatch->NextMiniBatch(miniSize, miniWords);

        // LOG(info)->info("Enqueuing batch of size {}", miniBatch->size());

        results.emplace_back(
          __GOD__.GetThreadPool().enqueue(
              [&__GOD__, miniBatch]{ return TranslationTask(__GOD__, miniBatch); }
              )
        );
      }

      maxiBatch.reset(new Sentences());
    }
  }

  // last batch
  if (maxiBatch->size()) {
    maxiBatch->SortByLength();
    while (maxiBatch->size()) {
      SentencesPtr miniBatch = maxiBatch->NextMiniBatch(miniSize, miniWords);
    //   LOG(info)->info("Enqueuing batch of size {}", miniBatch->size());
      results.emplace_back(
        __GOD__.GetThreadPool().enqueue(
            [&__GOD__, miniBatch]{ return TranslationTask(__GOD__, miniBatch); }
            )
      );
    }
  }

//   LOG(info)->info("Gathering histories.");
  // resort batch into line number order
  Histories allHistories;
  int historyIndex = 0;
  for (auto&& result : results) {
    // LOG(info)->info("Getting history {} ...", historyIndex);
    std::shared_ptr<Histories> histories = result.get();
    // LOG(info)->info("Got history {} ...", historyIndex);
    allHistories.Append(*histories);
    // LOG(info)->info("Saved history {}.", historyIndex);
  }
  allHistories.SortByLineNum();

//   LOG(info)->info("Saving to output.");
  // output
  for (size_t i = 0; i < allHistories.size(); ++i) {
    const History& history = *(allHistories.at(i).get());
    const Sentence& sentence = *(sentences_by_lineno[history.GetLineNum()]);
    saveHistoryToOutput(history,  sentence, __GOD__, _sentenceWords, _translations, _alignments);
  }
}

std::vector<int> getHardAlignment(Alignment_t& _alignment) {
    std::vector<int> _hardAlignment;
    for(size_t i = 0; i < _alignment.SourceSize; ++i) {
        float Max = 0.0;
        for(size_t j = 0; j < _alignment.TargetSize; ++j)
            if(_alignment.M[j * _alignment.SourceSize + i] > Max)
                Max = _alignment.M[j * _alignment.SourceSize + i];
        if(Max < 1e-5)
            continue;
        for(size_t j = 0; j < _alignment.TargetSize; ++j)
            _alignment.M[j * _alignment.SourceSize + i] /= Max;
    }
    _hardAlignment.resize(_alignment.TargetSize);
    for(size_t j = 0; j < _alignment.TargetSize; ++j) {
        float Max = _alignment.M[j * _alignment.SourceSize + 0];
        int ArgMax = 0;
        for(size_t i = 1; i < _alignment.SourceSize; ++i)
            if(_alignment.M[j * _alignment.SourceSize + i] > Max) {
                Max = _alignment.M[j * _alignment.SourceSize + i];
                ArgMax = i;
            }
        _hardAlignment[j] = ArgMax;
    }
    return _hardAlignment;
}

inline int minofVec(const std::vector<int>& _values)
{
    int Min = _values.at(0);
    for(size_t i = 1; i < _values.size(); ++i)
        if(Min > _values.at(i))
            Min = _values.at(i);
    return Min;
}

std::vector<std::tuple<int, int>> getCorrespondence(const std::string& _str1, const std::string& _str2)
{
#define DELETION     0
#define INSERTION    1
#define SUBSTITUTION 2

    const size_t L1 = _str1.size();
    const size_t L2 = _str2.size();
    const size_t S = L1 + 1;

    std::vector<int> D((L1 + 1) * (L2 + 1), 0);
    std::vector<int> Dij(3, 0);

    for(size_t i = 1; i <= L1; ++i)
        D[i] = 4 * i;
    for(size_t i = 1; i <= L2; ++i)
        D[i * S] = 4 * i;
    for(size_t i = 1; i <= L1; ++i) {
        for(size_t j = 1; j <= L2; ++j) {
            Dij[DELETION] = D[j * S + (i - 1)] + 1;
            Dij[INSERTION] = D[(j - 1) * S + i] + 1;
            Dij[SUBSTITUTION] = D[(j - 1) * S + (i - 1)];
            if(_str1.at(i - 1) != _str2.at(j - 1))
                Dij[SUBSTITUTION] += 1;
            D[j * S + i] = minofVec(Dij);
        }
    }

    std::vector<std::tuple<int, int>> Correspondence(_str2.size());

    size_t i = L1;
    size_t j = L2;
    size_t End = L1;
    while(i > 0 && j > 0) {
        Dij[DELETION] = D[j * S + (i - 1)] + 1;
        Dij[INSERTION] = D[(j - 1) * S + i] + 1;
        Dij[SUBSTITUTION] = D[(j - 1) * S + (i - 1)];
        if(_str1.at(i - 1) != _str2.at(j - 1))
            Dij[SUBSTITUTION] += 1;

        int Winner;
        if((_str1.at(i - 1) == ' ') != (_str2.at(j - 1) == ' ')) {
            if(_str1.at(i - 1) == ' ')
                Winner = DELETION;
            else
                Winner = INSERTION;
        } else {
            Winner = DELETION;
            if(Dij[Winner] > Dij[INSERTION])
                Winner = INSERTION;
            if(Dij[Winner] > Dij[SUBSTITUTION])
                Winner = SUBSTITUTION;
        }

        switch(Winner) {
        case DELETION:
            --i;
            break;
        case INSERTION:
            Correspondence[j - 1] = std::tuple<int, int>(-1, -1);
            --j;
            break;
        case SUBSTITUTION:
            Correspondence[j - 1] = std::tuple<int, int>(i - 1, End);
            End = i - 1;
            --i;
            --j;
            break;
        }
    }
    while(j > 0) {
        Correspondence[j - 1] = std::tuple<int, int>(-1, -1);
        --j;
    }

    return Correspondence;
}

void handleRequest(evhttp_request *_request, void *context)
{

    if (_request == nullptr)
        return;

    if (evhttp_request_get_command(_request) != EVHTTP_REQ_POST)
    {
        evhttp_send_error(_request, HTTP_BADMETHOD, "{ \"code\": 405, \"msg\": \"Method not allowed.\" }");
        return;
    }

    try {

        LOG(info)->info("Processing incoming request ...");

        boost::timer::cpu_timer timer;

        InputLines_t SrcAndOrigs = getTextLines(getPostBody(_request));

        std::vector<std::string>& Lines = std::get<0>(SrcAndOrigs);
        std::vector<std::string>& Origs = std::get<1>(SrcAndOrigs);

        God &__GOD__ = *(God *)context;

        std::vector<NBestTranslations_t> Translations;
        std::vector<std::vector<Alignment_t>> Alignments;
        std::vector<std::vector<std::string>> SourceTokens;
        
        LOG(info)->info("Found {} lines in {}", Lines.size(), timer.format());

        translate(Lines, __GOD__, SourceTokens, Translations, Alignments);

        LOG(info)->info("Translation took {}", timer.format());

        std::vector<std::vector<std::vector<std::string>>> AllPhrases;
        std::vector<std::vector<std::vector<int>>> AllAlignments;
        std::vector<std::vector<std::vector<std::tuple<int, int>>>> AllOrigAlignments;

        for(size_t Index = 0; Index < Lines.size(); ++Index) {
            std::vector<std::string>& Line = SourceTokens[Index];

            NBestTranslations_t& LineTranslations = Translations[Index];
            std::vector<Alignment_t>& LineAlignments = Alignments[Index];

            std::vector<std::string>& BestTranslation = LineTranslations[0];
            std::vector<int> WordMapping = getHardAlignment(LineAlignments[0]);
            // for(size_t iiii=0; iiii < LineAlignments[0].TargetSize; ++iiii) {
            //     for(size_t jjjj=0; jjjj < LineAlignments[0].SourceSize; ++jjjj)
            //         std::cerr << LineAlignments[0].M[iiii * LineAlignments[0].SourceSize + jjjj] << ", ";
            //     std::cerr << std::endl;
            // }
            // for(size_t iiii=0; iiii < BestTranslation.size(); ++iiii) 
            //     std::cerr << BestTranslation[iiii] << ": " << iiii << std::endl;
            // for(size_t iiii=0; iiii < WordMapping.size(); ++iiii) {
            //     int jjjj = WordMapping[iiii];
            //     std::cerr << BestTranslation[iiii] << ": " << iiii << " (" << Line[jjjj] << ")" << std::endl;
            // }

            std::vector<std::vector<std::string>> Phrases;
            std::vector<std::vector<int>> Alignments;
            std::vector<std::vector<std::tuple<int, int>>> OrigAlignments;
            size_t PhraseCount = 0;

            size_t LastSourceIndex = (size_t)-1;
            for(size_t TargetIndex = 0; TargetIndex < WordMapping.size(); ++TargetIndex) {
                std::string TargetWord = BestTranslation[TargetIndex];
                size_t SourceIndex = WordMapping[TargetIndex];
                if(TargetWord == UNK_STR) {
                    TargetWord = Line[SourceIndex];
                }
                if(SourceIndex == LastSourceIndex) {
                    Phrases[PhraseCount - 1][0] += std::string(" ") + TargetWord;
                } else {
                    Phrases.push_back({ TargetWord });
                    Alignments.push_back({ SourceIndex });
                    ++PhraseCount;
                }
                LastSourceIndex = SourceIndex;
            }

            for(size_t NBestIndex = 1; NBestIndex < LineTranslations.size(); ++NBestIndex) {
                std::vector<std::string>& TranslationCandidate = LineTranslations[NBestIndex];
                std::vector<int> WordMapping = getHardAlignment(LineAlignments[NBestIndex]);

                std::vector<std::string> CandidatesBySource;
                CandidatesBySource.resize(Line.size());
                for(size_t TargetIndex = 0; TargetIndex < WordMapping.size(); ++TargetIndex) {
                    auto TargetWord = TranslationCandidate[TargetIndex];
                    size_t SourceIndex = WordMapping[TargetIndex];
                    if(TargetWord == UNK_STR) {
                        TargetWord = Line[SourceIndex];
                    }
                    if(CandidatesBySource[SourceIndex].size() > 0)
                        CandidatesBySource[SourceIndex] += " ";
                    CandidatesBySource[SourceIndex] += TargetWord;
                }

                for(size_t PhraseIndex = 0; PhraseIndex < Phrases.size(); ++PhraseIndex) {
                    size_t SourceIndex = Alignments[PhraseIndex][0];
                    std::vector<std::string>& Candidates = Phrases[PhraseIndex];
                    bool Found = false;
                    for(size_t CandidateIndex = 0; CandidateIndex < Candidates.size(); ++CandidateIndex)
                        if(Candidates[CandidateIndex] == CandidatesBySource[SourceIndex]) {
                            Found = true;
                            break;
                        }
                    if(Found == false) {
                        Candidates.push_back(CandidatesBySource[SourceIndex]);
                    }
                }
            }

            AllPhrases.push_back(Phrases);
            AllAlignments.push_back(Alignments);
            AllOrigAlignments.push_back(OrigAlignments);
        }

        sendResults(_request, SourceTokens, AllPhrases, AllAlignments);

    } catch(std::exception& e) {
        std::cerr << "Exception during handling request: " << e.what() << std::endl;
        evhttp_send_error(_request, HTTP_INTERNAL, "{ \"code\": 500, \"msg\": \"Internal server error.\" }");
    }
}

typedef std::function<void(evhttp *)> SocketBinder_t;
void runEventLoop(std::atomic_bool &_done, SocketBinder_t _bind, God &__GOD__)
{
    std::unique_ptr<event_base, decltype(&event_base_free)> EventBase(
        event_base_new(), &event_base_free);
    if (!EventBase)
        throw std::runtime_error("Failed to create a new event base.");
    std::unique_ptr<evhttp, decltype(&evhttp_free)> EvHttp(
        evhttp_new(EventBase.get()),
        &evhttp_free);
    if (!EvHttp)
        throw std::runtime_error("Failed to create a new http event handler.");
    evhttp_set_gencb(EvHttp.get(), handleRequest, &__GOD__);
    _bind(EvHttp.get());
    while (_done == false)
    {
        event_base_loop(EventBase.get(), EVLOOP_NONBLOCK);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage: http_server <name> <address> <port> <path> [amun options, run with --help]" << std::endl;
        return 1;
    }

    SERVER_NAME = argv[1];
    const char *address = argv[2];
    uint16_t port = std::stoi(argv[3]);
    const char *path = argv[4];

    God __GOD__;
    __GOD__.Init(argc - 4, argv + 4);

    std::setvbuf(stdout, NULL, _IONBF, 0);
    std::setvbuf(stdin, NULL, _IONBF, 0);

    std::atomic_bool Done;

    auto ThreadDeleter = [&](std::thread *_thread) {
        Done = true;
        _thread->join();
        delete _thread;
    };

    typedef std::unique_ptr<std::thread, decltype(ThreadDeleter)> thread_ptr;
    typedef std::vector<thread_ptr> ThreadPool_t;

    std::exception_ptr InitException;
    evutil_socket_t Socket = -1;

    auto Binder = [&](evhttp *_evHttp) {
        if (Socket == -1)
        {
            evhttp_bound_socket *BoundSock = evhttp_bind_socket_with_handle(
                _evHttp,
                address,
                port);
            if (!BoundSock)
                throw std::runtime_error("Failed to bind server socket.");
            Socket = evhttp_bound_socket_get_fd(BoundSock);
            if (Socket == -1)
                throw std::runtime_error("Failed to get server socket descriptor for next instance.");
        }
        else
        {
            if (evhttp_accept_socket(_evHttp, Socket) == -1)
                throw std::runtime_error("Failed to bind server socket for new instance.");
        }
    };

    Done = false;

    ThreadPool_t Pool;
    unsigned int ThreadCount = std::thread::hardware_concurrency() - 1;
    for (int Index = 0; Index < ThreadCount; ++Index)
    {
        thread_ptr Thread(
            new std::thread([&]() {
                try
                {
                    runEventLoop(Done, Binder, __GOD__);
                }
                catch (...)
                {
                    InitException = std::current_exception();
                }
            }),
            ThreadDeleter);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        if (InitException != std::exception_ptr())
        {
            Done = true;
            std::rethrow_exception(InitException);
        }
        Pool.push_back(std::move(Thread));
    }

    /*std::cout << "Press `Enter` to quit." << std::endl;
    #std::cin.get();*/
    while(true) sleep(1);
    Done = true;
}
