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

## WECC240

The [WECC240 partition](../Partitions/WECC240/WECC240.partition.json) divides the case
into four regions. Each region compares its generator speeds against the
[monolithic simulation](WECC240.monolithic.solver.json) for a fault at bus 4004 (GRIZZLY).

From the same directory, run the monolithic case and then start one process per region.

```shell
$sim WECC240.monolithic.solver.json

$sim WECC240.SoCal.solver.json &
$sim WECC240.NorCal.solver.json &
$sim WECC240.Southwest.solver.json &
$sim WECC240.Northwest.solver.json
wait
```
