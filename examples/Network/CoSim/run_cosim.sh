#!/usr/bin/env sh

readonly _cosim_dir=$(dirname $(realpath $0))

cd ${_cosim_dir}

./CoSimServer -c 'TwoBusCoSimServer.case.json' -b 1 &
_srv_pid=$!

./CoSimClient -c 'TwoBusCoSimClient.case.json' -b 2
if [ $? -ne 0 ]; then
    _exit_code=1
fi

wait ${_srv_pid}
if [ $? -ne 0 ]; then
    _exit_code=1
fi

exit ${_exit_code}
