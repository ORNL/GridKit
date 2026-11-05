# Co-simulation

Both `DynamicSimulation` and `ContingencyAnalysis` are capable of exchanging
scalar signals with other GridKit instances over ZeroMQ. These examples
showcase this functionality in the following ways:

- Simple three-bus / two-bus network simulation.
- IEEE39 network bifurcation with one GridKit instance running most of
  the grid and another simulating a large generator and its dependents.

## Running the examples

These examples may be in the following way:

```sh
sh /path/to/build/examples/PhasorDynamics/CoSimulation/run_cosim.sh \
    --case ieee39
```

The available cases are `two-bus`, `three-bus`, `ieee39`, and
`ieee39-contingencies`. Ensure either the executables are on your `PATH`
or that you specify them with `--dynamic` and `--analysis`. By default,
the script uses a unique IPC endpoint.  Use the `--endpoint` flag to
override this.

## Building

Enable SUNDIALS, `GridKit_ENABLE_ZMQ=ON`. Build GridKit normally otherwise.
