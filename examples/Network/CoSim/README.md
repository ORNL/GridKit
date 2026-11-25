# Two-sided phasor dynamics co-simulation

This may be built with the following command:

```sh
nix develop
cmake -S . -B build-cosim -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DGridKit_ENABLE_SUNDIALS=ON -DGridKit_ENABLE_ZMQ=ON -DGridKit_ENABLE_ENZYME=ON
cmake --build build-cosim --target CoSimClient CoSimServer
cd build-cosim/examples/Network/CoSim
```

## Direct ZeroMQ exchange

```sh
bash run_cosim.sh three-bus
bash run_cosim.sh ieee39
bash run_cosim.sh ieee39 signals
```

The cilent and server exchange flat JSON objects over ZeroMQ. The client sends
bus voltage and the server returns accumulated bus current injections. Use
the same `--time-step` and `--final-time` on both peers when launching them
by hand. A `SignalVoltageOut` interface bus is used for the client and a
`SignalVoltageIn` interface bus is used on the server. Use `--output-bus`
to specify the interface bus for either the client or server.

## ARCADE usage

`CoSimClient` may alternately communicate over the
[ARCADE](https://github.com/sandialabs/ARCADE) semaphore protocol, with
`CoSimServer` communicating over ZeroMQ to the broker. To run ARCADE,
first copy the `IEEE39Arcade.input.json` to `input.json`. Then, run the
`DataBroker` executable. Once the broker reports that the external simulator
can start, run `bash run_cosim.sh ieee39 bus arcade` to use the broker
for communication between both sides of the simulation.
