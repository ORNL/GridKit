# **Generator Model**

The Generator model represents a dispatchable active- and reactive-power
injection with a quadratic active-power cost. It owns the dispatch variables,
adds their power contributions to the connected bus, and writes the accepted
dispatch and bus-current injection to the solution state.

Notes:
- Active and reactive generation are positive for injection into the connected
  bus.
- Dispatch, limits, and cost use per unit on the system power base. The machine
  power base `mva` is retained as model metadata and does not rescale these
  quantities.
- An offline generator remains in the model with both dispatch variables fixed
  to zero.

## Model Parameters

Symbol       | Units                 | Description                         | Typical Value | Note
-------------|-----------------------|-------------------------------------|---------------|------
$P_G^{\min}$ | [p.u.]               | Minimum active-power generation     |               | Optional; input name: `pmin`
$P_G^{\max}$ | [p.u.]               | Maximum active-power generation     |               | Optional; input name: `pmax`
$Q_G^{\min}$ | [p.u.]               | Minimum reactive-power generation   |               | Optional; input name: `qmin`
$Q_G^{\max}$ | [p.u.]               | Maximum reactive-power generation   |               | Optional; input name: `qmax`
$S^{\mathrm{base}}$ | [MVA]         | Generator machine power base        |               | Input name: `mva`
$c_0$        | [currency/hr]         | Constant generation-cost coefficient | 0             | Input name: `c0`
$c_1$        | [currency/(p.u. hr)]  | Linear generation-cost coefficient  | 0             | Input name: `c1`
$c_2$        | [currency/(p.u.^2 hr)] | Quadratic generation-cost coefficient | 0             | Input name: `c2`

### Parameter Validation

Invalid Generator parameter data are rejected by the following checks:

```math
\begin{aligned}
  &S^{\mathrm{base}} \in \mathbb{R}\ \text{and finite},
    \qquad S^{\mathrm{base}} > 0 \\
  &P_G^{\min}, P_G^{\max}, Q_G^{\min}, Q_G^{\max}
    \in \mathbb{R}\ \text{and finite when supplied} \\
  &P_G^{\min} \le P_G^{\max}
    \quad\text{when both active-power limits are supplied} \\
  &Q_G^{\min} \le Q_G^{\max}
    \quad\text{when both reactive-power limits are supplied} \\
  &c_0, c_1, c_2 \in \mathbb{R}\ \text{and finite}
\end{aligned}
```

## Model Inputs

Input    | Symbol | Units  | Description                         | Default
---------|--------|--------|-------------------------------------|--------
`online` | $u$    | [-]    | Generator connection status         | 1

The `online` input is constant while a solve is running.

### Model Derived Parameters

Define the connection indicator as:

```math
u =
\begin{cases}
  1, & \text{when online is true or omitted} \\
  0, & \text{when online is false}
\end{cases}
```

For an online generator, missing limits are treated as unbounded. For an
offline generator, the effective lower and upper bounds are both zero:

```math
\begin{aligned}
  \underline{P}_G &=
  \begin{cases}
    P_G^{\min}, & u=1 \text{ and } P_G^{\min} \text{ is supplied} \\
    -\infty,    & u=1 \text{ and } P_G^{\min} \text{ is omitted} \\
    0,          & u=0
  \end{cases} \\
  \overline{P}_G &=
  \begin{cases}
    P_G^{\max}, & u=1 \text{ and } P_G^{\max} \text{ is supplied} \\
    \infty,     & u=1 \text{ and } P_G^{\max} \text{ is omitted} \\
    0,          & u=0
  \end{cases} \\
  \underline{Q}_G &=
  \begin{cases}
    Q_G^{\min}, & u=1 \text{ and } Q_G^{\min} \text{ is supplied} \\
    -\infty,    & u=1 \text{ and } Q_G^{\min} \text{ is omitted} \\
    0,          & u=0
  \end{cases} \\
  \overline{Q}_G &=
  \begin{cases}
    Q_G^{\max}, & u=1 \text{ and } Q_G^{\max} \text{ is supplied} \\
    \infty,     & u=1 \text{ and } Q_G^{\max} \text{ is omitted} \\
    0,          & u=0
  \end{cases}
\end{aligned}
```

## Model Variables

### Internal Variables

Symbol | Units  | Description                       | Note
-------|--------|-----------------------------------|------
$P_G$  | [p.u.] | Active-power generation           | Enum name: `PG`; bounded by $\underline{P}_G$ and $\overline{P}_G$
$Q_G$  | [p.u.] | Reactive-power generation         | Enum name: `QG`; bounded by $\underline{Q}_G$ and $\overline{Q}_G$

### External Variables

Symbol | Units  | Description                    | Note
-------|--------|--------------------------------|------
$V_M$  | [p.u.] | Connected-bus voltage magnitude | Enum name: `VM`; owned by connected bus
$V_A$  | [rad]  | Connected-bus voltage angle     | Enum name: `VA`; owned by connected bus

## Wiring

Port  | Type | Description
------|------|------------
`bus` | Bus  | Connected bus that owns terminal voltage variables and power-balance constraints

### Configuration Validation

Each generator `id` must be nonempty and unique across all devices. Its `bus`
value must identify an existing bus.

## Model Constraints

### Internal Constraints

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

### External Constraints

Symbol          | Units  | Description                           | Note
----------------|--------|---------------------------------------|------
$\mathrm{DIVP}$ | [p.u.] | Connected-bus active-power balance    | Enum name: `DIVP`; owned by connected bus
$\mathrm{DIVQ}$ | [p.u.] | Connected-bus reactive-power balance  | Enum name: `DIVQ`; owned by connected bus

The generator adds its solved dispatch to the connected-bus balances:

```math
\begin{aligned}
  \mathrm{DIVP} &\mathrel{+}= P_G \\
  \mathrm{DIVQ} &\mathrel{+}= Q_G
\end{aligned}
```

## Objective Contribution

The online generator contributes its quadratic active-power cost. An offline
generator contributes zero cost, including no constant term:

```math
f_G = u\left(c_0 + c_1 P_G + c_2 P_G^2\right)
```

The corresponding objective-gradient contribution is:

```math
\begin{aligned}
  \dfrac{\partial f_G}{\partial P_G}
    &= u\left(c_1 + 2c_2P_G\right) \\
  \dfrac{\partial f_G}{\partial Q_G}
    &= 0
\end{aligned}
```

## Initialization

Input | Symbol    | Units  | Description                       | Default
------|-----------|--------|-----------------------------------|--------
`p`   | $P_{G,0}$ | [p.u.] | Initial active-power generation   | 0
`q`   | $Q_{G,0}$ | [p.u.] | Initial reactive-power generation | 0

The companion generator state supplies the initialization inputs above. Any
supplied values define only the initial point; the accepted values are
determined by the optimization.

Any supplied initial dispatch must satisfy:

```math
P_{G,0}, Q_{G,0} \in \mathbb{R}\ \text{and finite}
```

The dispatch variables initialize from the companion state when the generator
is online and initialize to zero otherwise:

```math
\begin{aligned}
  P_G &= u P_{G,0} \\
  Q_G &= u Q_{G,0}
\end{aligned}
```

Missing `p` or `q` values and a missing generator state record use zero. The
initialized values are not clamped to the adjustable limits because they form
an optimization starting point.

## Model Outputs

Output | Units  | Description                                  | Note
-------|--------|----------------------------------------------|------
`p`    | [p.u.] | Accepted active-power generation             | Written to the generator device state
`q`    | [p.u.] | Accepted reactive-power generation           | Written to the generator device state
`ir`   | [p.u.] | Bus-current injection, real component        | Written under the generator `id` at the connected bus
`ii`   | [p.u.] | Bus-current injection, imaginary component   | Written under the generator `id` at the connected bus

For the accepted connected-bus voltage, define:

```math
\begin{aligned}
  V_r &= V_M \cos V_A \\
  V_i &= V_M \sin V_A
\end{aligned}
```

The generator current injected into the bus follows $S=VI^{\ast}$:

```math
\begin{aligned}
  I_r &= \dfrac{P_G V_r + Q_G V_i}{V_M^2} \\
  I_i &= \dfrac{P_G V_i - Q_G V_r}{V_M^2}
\end{aligned}
```

For an online generator, state output requires $V_M>0$. An offline generator
writes zero current without evaluating the division above. The generator
device record is created when it is missing. Existing status fields and
unrelated state data are preserved.
