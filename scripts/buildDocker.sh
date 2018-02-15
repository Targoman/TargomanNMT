#!/bin/bash

pushd $(dirname $0) >/dev/null
cd ..

docker tag targoman/amun_rest_server:latest targoman/amun_rest_server:$(date +"%y%m%d_%H%M%S")
docker build -t targoman/amun_rest_server -f scripts/docker/Dockerfile.server .

popd > /dev/null
