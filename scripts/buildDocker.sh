#!/bin/bash

pushd $(dirname $0) >/dev/null
cd ..

docker rmi --force targoman/amun_rest_server
docker build -t targoman/amun_rest_server -f scripts/docker/Dockerfile.server .

popd > /dev/null
