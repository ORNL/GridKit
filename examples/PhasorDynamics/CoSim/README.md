# CoSimulation

`DynamicSimulation` simulates two partitions at one interfacing `Bus` at a shared `dt_monitor`
when both solver files have a `connect` event for that bus.

The IEEE39 case is split at bus 39:

| Partition | Contents | Interface bus |
| --- | --- | --- |
| A | Network and generators 30 to 38 | `SignalVoltageIn` |
| B | Bus 39 generator, exciter and stabilizer | `SignalVoltageOut` |

The reference files are the generator speeds of the monolithic
[IEEE39 validation](../Validation/IEEE39/README.md) study.

## Tutorial

Start one process per partition, in the same directory, within 10 seconds of each other:
```shell
cd build/examples/PhasorDynamics/CoSim
sim=../../../application/PhasorDynamics/DynamicSimulation

$sim IEEE39.A.solver.json &
$sim IEEE39.B.solver.json
```

Each process writes its own output and error against the reference obtained from the monolithic solution.
