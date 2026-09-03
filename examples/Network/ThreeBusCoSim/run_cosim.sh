#!/usr/bin/env sh

readonly _cosim_dir=$(dirname $(realpath $0))

cd ${_cosim_dir}

./CoSimServer &

./CoSimClient

wait
