# Distributed-line geometry

`overhead.line.json` and `north-american.catalog.json` come from
[GridWorkbench](https://github.com/abirchfield/GridWorkbench), commit
`3362db50c065c79b373f50029a03eb15fd9aa9b6`, under
`examples/emt/geometry/345kv-horizontal.line.json` and
`examples/emt/geometry/north-american.catalog.json`.
The source describes these conductor, tension, and earth data as illustrative.

The study script uses the supplied twin Drake ACSR bundles, shield wires,
catenary geometry, and 0.01 S/m earth conductivity. It applies cyclic phase
transposition to GridWorkbench's reduced phase matrices before fitting.
The line's configured path length is overridden by each study edge length.

See [the study instructions](../../../examples/EMT/Distributed/README.md)
for fitting, simulation, and plotting. Generated coefficients and results
remain local; these files preserve the input geometry and provenance.
