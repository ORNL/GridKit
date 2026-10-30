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

## Gas Laws

A GasBranch takes the gas law as $a$, $b_1$, and $b_2$, with density
$\rho(p) = p\,(b_1 + b_2 p)/a^2$. GasTranSim derives them from the gas
temperature $T$ [K] and specific gravity $G$ [-]:

```math
a = \sqrt{\dfrac{R_u T}{M_\mathrm{air} G}}
```

with $R_u = 8.314$ J/(mol·K) and $M_\mathrm{air} = 0.02896$ kg/mol.

Gas law     | $b_1$                   | $b_2$ [1/Pa]
------------|-------------------------|--------------
Ideal       | 1                       | 0
Simple CNGA | 1.00300865              | 2.96848838e-8
Full CNGA   | $1 + c\,p_\mathrm{atm}$ | $c$

Here $c = 344400 \cdot 10^{1.785G} / \left(6894.75729\,(1.8T)^{3.825}\right)$
and $p_\mathrm{atm} = 101350$ Pa.

## References

- V. Gyrya and A. Zlotnik, [“An explicit staggered-grid method for numerical
  simulation of large-scale natural gas pipeline
  networks,”](https://doi.org/10.1016/j.apm.2018.07.051) *Applied Mathematical
  Modelling*, 2019.
- [GasTranSim.jl](https://github.com/kaarthiksundar/GasTranSim.jl)
