# Fitting Applications

Applications for rational approximation of sampled frequency responses,
driven by the
[Optimization solver family](../../GridKit/Solver/Optimization/VectorFitting/README.md).

`VectorFitting` is the model-agnostic driver: the solver specification names
a sampled-response CSV and the fitting options; the application writes the
fitted rational model JSON, an optional state-space realization JSON, and
prints the fit statistics.

EMT line fitting, which produces its samples in-process and adds delay
extraction, is the separate `UniversalLineModel` application in
[`application/EMT`](../EMT).
