# **Bus Model**

The Bus model represents a network connection point in the optimal power-flow
formulation. Each bus owns voltage-magnitude and voltage-angle variables and
the active- and reactive-power balance constraints. Connected components add
their signed power-injection contributions directly to these constraints.

Notes:
- Power entering the bus has positive sign, and power leaving the bus has
  negative sign.
- Each bus is identified by its `id` and `number`.
- A bus with class `Slack` provides the voltage-angle reference.

## Model Parameters

Symbol            | Units  | Description                         | Typical Value | Note
------------------|--------|-------------------------------------|---------------|------
$n$               | [-]    | Bus number                          |               | Input name: `number`
$V_\mathrm{base}$   | [kV] | Nominal bus voltage                 |               | Input name: `kv`
$V_M^{\min}$      | [p.u.] | Minimum voltage magnitude           |               | Optional; input name: `vmin`
$V_M^{\max}$      | [p.u.] | Maximum voltage magnitude           |               | Optional; input name: `vmax`

### Parameter Validation

Invalid Bus parameter data are rejected by the following checks:

```math
\begin{aligned}
  &V_\mathrm{base} \in \mathbb{R}\ \text{and finite},
    \qquad V_\mathrm{base} > 0 \\
  &V_M^{\min}, V_M^{\max} \in \mathbb{R}\ \text{and finite}
    \quad\text{when supplied} \\
  &V_M^{\min} \ge 0
    \quad\text{when } V_M^{\min} \text{ is supplied} \\
  &V_M^{\max} > 0
    \quad\text{when } V_M^{\max} \text{ is supplied} \\
  &V_M^{\min} \le V_M^{\max}
    \quad\text{when both limits are supplied}
\end{aligned}
```

### Configuration Validation

Bus `id` values must be nonempty and unique. Bus `number` values must be
unique, with $n\in\mathbb{Z}_{\ge 0}$. The connected system must contain
exactly one `Slack` bus, and the closed-branch topology must connect every bus
to that `Slack` bus.

## Model Variables

### Internal Variables

Symbol | Units  | Description       | Note
-------|--------|-------------------|------
$V_M$  | [p.u.] | Voltage magnitude | Enum name: `VM`
$V_A$  | [rad]  | Voltage angle     | Enum name: `VA`

Missing voltage-magnitude limits are treated as unbounded. For a standard
`Bus`, the variable bounds are:

```math
\begin{aligned}
  V_M^{\min} &\le V_M \le V_M^{\max} \\
  -\infty &\le V_A \le \infty
\end{aligned}
```

When `vmin` or `vmax` is omitted, the corresponding bound above is
$-\infty$ or $\infty$, respectively. For a `Slack` bus, the lower and upper
voltage-angle bounds are equal.

### External Variables

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

## Model Constraints

### Internal Constraints

Symbol             | Units  | Description                    | Note
-------------------|--------|--------------------------------|------
$\mathrm{DIVP}_i$ | [p.u.] | Active-power balance at bus $i$ | Enum name: `DIVP`; equality constrained to zero
$\mathrm{DIVQ}_i$ | [p.u.] | Reactive-power balance at bus $i$ | Enum name: `DIVQ`; equality constrained to zero

Let $\mathcal{D}_i$ be the set of device terminals connected to bus $i$.
Each connected component contributes power oriented entering the bus:

```math
\begin{aligned}
  0 &= \mathrm{DIVP}_i
     := \sum_{d \in \mathcal{D}_i} P_{d \rightarrow i} \\
  0 &= \mathrm{DIVQ}_i
     := \sum_{d \in \mathcal{D}_i} Q_{d \rightarrow i}
\end{aligned}
```

The system resets both balance values to zero before connected components add
their contributions.

### External Constraints

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

## Initialization

Let $V_{r,0}$ and $V_{i,0}$ be the real and imaginary voltage components from
the required companion state record named `bus_id_<number>`. Initialization
requires:

```math
\begin{aligned}
  &V_{r,0}, V_{i,0} \in \mathbb{R}\ \text{and finite} \\
  &V_{r,0}^2 + V_{i,0}^2 > 0
\end{aligned}
```

The initial polar voltage is:

```math
\begin{aligned}
  V_{M,0} &= \sqrt{V_{r,0}^2 + V_{i,0}^2} \\
  V_{A,0} &= \operatorname{atan2}\!\left(V_{i,0}, V_{r,0}\right)
\end{aligned}
```

For a `Slack` bus, the initialized angle defines the fixed reference bounds:

```math
\underline{V}_A = \overline{V}_A = V_{A,0}
```

The initial voltage magnitude is not clamped to its limits. An initial point
outside the adjustable bounds remains a valid, infeasible starting point.

## Model Outputs

Output | Units  | Description                         | Note
-------|--------|-------------------------------------|------
`vr`   | [p.u.] | Solved bus voltage, real component  | Written to `bus_id_<number>`
`vi`   | [p.u.] | Solved bus voltage, imaginary component | Written to `bus_id_<number>`

The solved Cartesian voltage is:

```math
\begin{aligned}
  V_r &= V_M \cos V_A \\
  V_i &= V_M \sin V_A
\end{aligned}
```

Existing bus-injection records and unrelated state data are preserved.
