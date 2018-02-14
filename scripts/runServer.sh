#!/bin/bash

ServerName=$1
SourceLang=$2
DestLang=$3
GPU=$4
Model=$5

if [ -z "$ServerName" -o -z "$SourceLang" -o -z "$DestLang" -o -z "$GPU" -o -z "$Model" ];then
    echo "USAGE: $0 SourceLang DestLang GPU ModelPath"
    exit 1;
fi

./amun_rest_server "$ServerName" 0.0.0.0 5000 /translate \
    -m "/nmt/models/$Model/model.npz" \
    -s "/nmt/models/$Model/vocab.$SourceLang.yml" \
    -t "/nmt/models/$Model/vocab.$DestLang.yml" \
    -d $GPU \
    -b 12 --n-best -n --mini-batch 10 --maxi-batch 1000 \
    -u --return-alignment
