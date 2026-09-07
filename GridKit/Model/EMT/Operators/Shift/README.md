# Shift Operators

Shift operators represent signal propagation in time.

`Shift` supplies the common input/output signal, frequency-response, and
constant or harmonic prehistory initialization interface. Implementations
participate in the Component allocation, sparse Jacobian, accepted-history,
and discontinuity lifecycle. `Delay` owns its delayed algebraic outputs;
`Propagation` owns its rational and delay submodels and exposes their summed
output as a computed signal.

## Models

- [Delay](Delay/README.md)
- [Propagation](Propagation/README.md)
