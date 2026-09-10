# GENSAL Two-Bus Validation

This validation case compares GridKit's GENSAL response against an external
PowerWorld reference on a two-bus system.

The case contains a GENSAL machine on Bus 1, a 100 MW resistive load on Bus 2,
a lossless tie line with X = 0.1 p.u., and a temporary bus fault on Bus 2. The
monitored machine states are compared against `GENSAL.ref.csv`.

## Trajectory Comparison

![GENSAL validation trajectory](Gensal_validation.png)

## Error

![GENSAL validation error](Gensal_validation_error.png)
