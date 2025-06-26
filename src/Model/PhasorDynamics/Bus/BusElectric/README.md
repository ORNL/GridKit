# Electric Bus Model

A bus is a point of interconnection of electrical devices. Bus component model
also plays a key role in coupling system components. Each bus $`k`$ owns two
variables -- real and imaginary voltage $`V_{rk}`$ and $`V_{ik}`$,
respectively. The bus also owns current balance residual equations for real
and imaginary currents coming into the bus $`I_{rk}`$ and $`I_{ik}`$,
respectively. While th bus model owns current residuals, it _does not compute_
them. Instead, each component connected to the bus is adding its contribution
to the residual. The bus will merely initialize the residual to zero each time
numerical integrator requests residual evaluation.

**Sign Convention**

Current entering the bus has positive and current exiting the bus negative
sign.

<div align="center">
   <img align="center" src="../../../../docs/Figures/bus_variables.jpg">
   
  Figure 1: Needs to be changed to represent current balance instead of power
  balance.
</div>

**Other Parameters**
Buses are uniquely defined by their ID (number or name). Besides, each bus
should have associated Nominal Voltage value.

## Electric Bus Types

- Network Bus  (See [BusNetwork](./BusNetwork/))
- Infinite Bus (See [BusInfinite](./BusInfinite/))