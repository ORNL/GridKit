# Partitioned ACTIVSg70k

This example runs the
[ACTIVSg70k case](../../../../cases/PhasorDynamics/ACTIVSg70k/README.md) for
10 seconds with a bus fault applied at 1 s and cleared at 1.15 s, using the
IDA settings of the monolithic DOE study. The case is split into 64 regions
advanced by multicolor Gauss-Seidel on 8 OpenMP threads; build with
`GridKit_ENABLE_OPENMP=ON` and set `partition.threads` to the number of
performance cores.

Remove `partition` to run the same study monolithically, and add
`output_file` to record generator speeds.
