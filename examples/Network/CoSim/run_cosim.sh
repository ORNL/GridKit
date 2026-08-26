#!/bin/bash

./CoSimServer -c 'TwoBusCoSimServer.case.json' -b 1 &

./CoSimClient -c 'TwoBusCoSimClient.case.json' -b 2

wait
