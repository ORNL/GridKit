# ContingencyAnalysis

`ContingencyAnalysis` runs a PhasorDynamics study for each configured bus fault.
When OpenMP is available, the fault studies run through the OpenMP path; otherwise
the application uses asynchronous tasks.

## Inputs

- [Application Input Format](../../generated/application-input-format.md)
- [System Model Input Format](../../models/generated/phasor-dynamics/input-format.md)

## Source

- `application/PhasorDynamics/ContingencyAnalysis.cpp`
