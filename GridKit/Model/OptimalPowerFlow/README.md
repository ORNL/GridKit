# Optimal Power Flow

AC optimal power flow in Cartesian coordinates. It finds the least-cost
operating point within the component limits and writes it as a
[state](../STATE.md) that PhasorDynamics starts from.

## Formulation

```math
\begin{aligned}
\min_x \quad & \sum_g a_g \left(c_{0,g} + c_{1,g} P_g + c_{2,g} P_g^2\right) \\
\text{s.t.} \quad & \sum_{c \to k} P_{c,k}(x) = 0, \quad \sum_{c \to k} Q_{c,k}(x) = 0 && \text{each bus } k \text{ that is not infinite} \\
& P_{\ell,k}^2 + Q_{\ell,k}^2 \le \left(S_\ell^{\max}\right)^2 && \text{each terminal } k \text{ of a rated branch } \ell \\
& \left(V_k^{\min}\right)^2 \le V_{\mathrm{r}k}^2 + V_{\mathrm{i}k}^2 \le \left(V_k^{\max}\right)^2 \\
& P_g^{\min} \le P_g \le P_g^{\max}, \quad Q_g^{\min} \le Q_g \le Q_g^{\max} \\
& V_{\mathrm{r}}^\mathrm{state} V_{\mathrm{i}} - V_{\mathrm{i}}^\mathrm{state} V_{\mathrm{r}} = 0 && \text{reference bus}
\end{aligned}
```

- $P_{c,k}$ and $Q_{c,k}$ are the power component $c$ injects into bus $k$.
  Every sign in this family is positive into the bus.
- $a_g$ is one for an online generator and zero for an offline one.
- The reference bus is the lowest-numbered bus. With an infinite bus there is
  no reference bus. An infinite bus keeps its state voltage, has no balance
  rows, and supplies any mismatch.
- A constraint without a finite bound gets no row.

## Components

Component   | Variables            | Constraints                      | Kernels
------------|----------------------|----------------------------------|----------
`Bus`       | $V_\mathrm{r}$, $V_\mathrm{i}$ | Power balance, $\lvert V \rvert^2$, and angle reference | Quadratic
`Branch`    | Terminal $V_\mathrm{r}$, $V_\mathrm{i}$ | $\lvert S_1 \rvert^2$, $\lvert S_2 \rvert^2$ and power into both buses | Quartic
`Generator` | $P$, $Q$ and terminal $V_\mathrm{r}$, $V_\mathrm{i}$ | Power into the bus  | Quadratic cost
`Load`      | Terminal $V_\mathrm{r}$, $V_\mathrm{i}$ | Power into the bus, fixed by the state | Constant
`Shunt`     | Terminal $V_\mathrm{r}$, $V_\mathrm{i}$ | Power into the bus               | Quadratic

Local variables are the internal variables followed by $V_\mathrm{r}$ and
$V_\mathrm{i}$ of each terminal bus. Local constraints are the internal
constraints followed by $P$ and $Q$ into each terminal bus. `SystemModel` maps both to global indices, so a
device adds its injection to the balance rows of its buses.

## Derivatives

Each model supplies two kernels over its local variables, `objective(x)` and
`constraints(x, g)`. `ComponentModel` differentiates them with Enzyme:

Driver   | Derivative                                           | Mode
---------|------------------------------------------------------|---------------------
`DfDx`   | Objective gradient                                   | Reverse
`DgDx`   | Constraint Jacobian                                  | Forward, one sweep per local variable
`D2LDx2` | Lower triangle of the Hessian of $\sigma f + \lambda^T g$ | Forward over reverse

The Jacobian and Hessian stores are sparse. Enzyme's sparsity analysis stores
every structural entry, numerical zeros included, so the patterns set in
`allocate()` hold at every point and for every multiplier. `Load` has a
constant kernel and no entries.

The `DependencyTracking::Variable` instantiations give the same Jacobian and
objective gradient from dependencies. Tests compare the two.

## Data

`SystemModelData` holds the network. The `OptimalDispatch` application builds
it from a PhasorDynamics case, and `applyMatpowerData` adds limits and costs
from a [MATPOWER case](https://matpower.org/docs/ref/matpower/caseformat.html)
of the same network:

- Buses match by number and take `VMIN` and `VMAX`.
- The in-service branches between two buses match the branches between them
  in order and take a nonzero `RATE_A`.
- The in-service generators at a bus match the generators at that bus in
  order and take `PMIN`, `PMAX`, `QMIN`, `QMAX`, and their polynomial cost of
  at most second order. MATPOWER generators at a bus without generators are
  static injections that the loads already carry.
- Values in MW, Mvar, and MVA and costs in P [MW] become system base values.

`parseMatpowerData` reads every numeric `mpc.<name>` matrix. Rows end at `;` or
at the end of a line, so both MATPOWER and PowerWorld exports read.

## State

`SystemModel` reads the starting point, the load demand, and the device settings
`online`, `open`, `tap`, and `phase` from a state. `solutionState()` returns
that state with the solved bus voltages and the terminal currents of every
device, computed from the power into each terminal bus.
