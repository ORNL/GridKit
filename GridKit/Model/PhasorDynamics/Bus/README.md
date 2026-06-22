# Bus Model

A bus is a point of interconnection of electrical devices. The bus component
model also plays a key role in coupling system components. Each bus $k$ owns
two variables: real voltage and imaginary voltage, denoted as $V_{rk}$ and
$V_{ik}$, respectively. The bus also owns current-balance residual equations for
real and imaginary currents entering the bus, denoted as $I_{rk}$ and $I_{ik}$,
respectively. While the bus model owns current residuals, it _does not compute_
them. Instead, each component connected to the bus adds its contribution to the
residual. The bus initializes the residual to zero each time the numerical
integrator requests residual evaluation.

## Sign Convention

Current entering the bus has positive sign, and current exiting the bus has
negative sign.

![](../../../../docs/Figures/bus_variables.jpg)

Figure 1: Needs to be changed to represent current balance instead of power

balance.
**Other Parameters**
Buses are uniquely defined by their ID (number or name). Besides, each bus
should have associated Nominal Voltage value.
