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
./run_cosim.sh ieee39-fault      # IEEE39 bus 16 fault, 10-second simulation

# test case formulations
ctest -R 'CoSim.*(Example|Test)$' --output-on-failure
```

The fault case uses the IEEE39 validation disturbance: R = 0 and X = 0.01 pu
at bus 16, applied at 1.0 s and cleared at 1.15 s. The grid partition applies
the fault. Outputs are `IEEE39Grid.csv` and `IEEE39Generator.csv`.
