# Co-simulation example

`CoSimClient` and `CoSimServer` exchange flat JSON objects over ZeroMQ
REQ/REP. The client sends bus voltage and the server returns accumulated
bus current injections.

The supplied interface buses are `SignalVoltageOut` on the client and
`SignalVoltageIn` on the server. Their input ports must be connected to
the signals in the case file. Output signal ports are optional when reading
the bus directly.

The co-simulated cases may be run in the following ways:

```sh
./run_cosim.sh                   # three-bus using bus outputs
./run_cosim.sh ieee39            # IEEE39 using bus outputs
./run_cosim.sh ieee39 signals    # IEEE39 using signal outputs

# test case formulations
ctest -R 'CoSim.*(Example|Test)$' --output-on-failure
```

