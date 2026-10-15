# EconomicDispatch

These examples solve the AC optimal power flow of MATPOWER cases with Ipopt,
using the exact sparse derivatives of the `OptimalPowerFlow` model.

| Case | Buses | Optimal cost [\$/h] |
| --- | --- | --- |
| [IEEE39](../../../cases/OptimalPowerFlow/IEEE39/README.md) | 39 | 41864.18 |
| [Hawaii](../../../cases/OptimalPowerFlow/Hawaii/README.md) | 37 | 4151.35 |
| [ACTIVSg2000](../../../cases/OptimalPowerFlow/ACTIVSg2000/README.md) | 2000 | 1228892.07 |

The tests check each cost against PYPOWER 5.1.21, a port of MATPOWER, to a
relative 1e-8.
