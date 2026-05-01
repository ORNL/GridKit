# Bus Model

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

## Actions

This component accepts the following runtime events via `apply(Action)`. See
[EVENTS.md](../EVENTS.md) for the dispatch model and JSON schema.

JSON  | Params   | Effect
----------------|----------|--------------------------------------------------------------------
`"fault"`       | `R`, `X` | Stores the fault impedance and engages the fault ($U = 1$).
`"clear"`       | none     | Disengages the fault ($U = 0$). Stored impedance is unchanged.

While the fault is engaged, the bus current balance gains the contributions:

``` math
\begin{aligned}
  G_f &= \dfrac{R}{R^2 + X^2}, \quad B_f = -\dfrac{X}{R^2 + X^2} \\
  \Delta I_{rk} &= U(-G_f V_{rk} + B_f V_{ik}) \\
  \Delta I_{ik} &= U(-B_f V_{rk} - G_f V_{ik})
\end{aligned}
```

The `percent` parameter of `Fault` is ignored by `Bus`. The fault gate $U$
is implemented as a mask (0 or 1) on the fault residual term,
so the Jacobian sparsity pattern is fixed across `Fault`/`Clear` cycles.
