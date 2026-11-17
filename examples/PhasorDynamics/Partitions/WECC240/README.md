# WECC240 Partitioned

This example splits the reusable
[WECC240 case](../../../../cases/PhasorDynamics/WECC240/README.md)
into four subsystems with `WECC240.partition.json` and runs the same fault
study as the [validation example](../../Validation/WECC240/README.md).

In this study the partition only organizes the model into subsystems, and a
single IDA integrator still solves the whole system.

The study has no reference file. Monitor columns of a partitioned system are in
a different order, and the reference comparison goes by column position.

When SUNDIALS provides the partitioned DAE integrator,
`WECC240.integrator.solver.json` runs two steps without events through
`PartitionedDynamicSimulation`, with one partition per subsystem. The study only
checks that the run completes.
