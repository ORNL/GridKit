# Gas Network Models

Gas network models represent isothermal pipeline flow. Pressure plays the role
of bus voltage, and mass flow plays the role of current.

## Conventions

- Quantities are in SI units.
- Mass flow entering a bus is positive.
- A synchronous machine draws its fuel from a gas bus through its `GasInputPm`
  port.

## Types

- [GasBus](GasBus/README.md): Junction mass balance
- [GasBusInfinite](GasBusInfinite/README.md): Fixed-pressure supply
- [GasBranch](GasBranch/README.md): Pipe

## References

- V. Gyrya and A. Zlotnik, [“An explicit staggered-grid method for numerical
  simulation of large-scale natural gas pipeline
  networks,”](https://doi.org/10.1016/j.apm.2018.07.051) *Applied Mathematical
  Modelling*, 2019.
- [GasTranSim.jl](https://github.com/kaarthiksundar/GasTranSim.jl)
