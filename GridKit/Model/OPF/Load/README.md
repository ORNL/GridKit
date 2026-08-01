# **Load Model**

The Load model represents a fixed active- and reactive-power injection read
from the companion state. It owns no optimization variables or constraints and
adds its signed power contribution directly to the connected-bus balances.

Notes:
- Load demand uses the network injection convention. Consumed active and
  reactive power are stored as negative `p` and `q` values.
- The base formulation does not dispatch load. The optional limits validate
  the fixed operating state and do not become optimization-variable bounds.
- An offline load contributes zero power and current without changing the
  stored `p` and `q` values.

## Model Parameters

Symbol         | Units  | Description                       | Typical Value | Note
---------------|--------|-----------------------------------|---------------|------
$P_L^{\min}$  | [p.u.] | Minimum signed active-power injection |               | Optional; input name: `pmin`
$P_L^{\max}$  | [p.u.] | Maximum signed active-power injection |               | Optional; input name: `pmax`
$Q_L^{\min}$  | [p.u.] | Minimum signed reactive-power injection |             | Optional; input name: `qmin`
$Q_L^{\max}$  | [p.u.] | Maximum signed reactive-power injection |             | Optional; input name: `qmax`

### Parameter Validation

Invalid Load parameter data are rejected by the following checks:

```math
\begin{aligned}
  &P_L^{\min}, P_L^{\max}, Q_L^{\min}, Q_L^{\max}
    \in \mathbb{R}\ \text{and finite when supplied} \\
  &P_L^{\min} \le P_L^{\max}
    \quad\text{when both active-power limits are supplied} \\
  &Q_L^{\min} \le Q_L^{\max}
    \quad\text{when both reactive-power limits are supplied}
\end{aligned}
```

## Model Inputs

Input    | Symbol | Units  | Description                           | Default
---------|--------|--------|---------------------------------------|--------
`online` | $u$    | [-]    | Load connection status                | 1
`p`      | $p$    | [p.u.] | Fixed signed active-power injection   | Required
`q`      | $q$    | [p.u.] | Fixed signed reactive-power injection | Required

These inputs are constant while a solve is running. The companion device state
must supply `p` and `q`, and the fixed operating inputs must satisfy:

```math
\begin{aligned}
  &p, q \in \mathbb{R}\ \text{and finite} \\
  &P_L^{\min} \le p
    \quad\text{when } P_L^{\min} \text{ is supplied} \\
  &p \le P_L^{\max}
    \quad\text{when } P_L^{\max} \text{ is supplied} \\
  &Q_L^{\min} \le q
    \quad\text{when } Q_L^{\min} \text{ is supplied} \\
  &q \le Q_L^{\max}
    \quad\text{when } Q_L^{\max} \text{ is supplied}
\end{aligned}
```

### Model Derived Parameters

Define the connection indicator as:

```math
u =
\begin{cases}
  1, & \text{when online is true or omitted} \\
  0, & \text{when online is false}
\end{cases}
```

## Model Variables

### Internal Variables

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

### External Variables

Symbol | Units  | Description                       | Note
-------|--------|-----------------------------------|------
$V_M$  | [p.u.] | Connected-bus voltage magnitude   | Enum name: `VM`; owned by connected bus
$V_A$  | [rad]  | Connected-bus voltage angle       | Enum name: `VA`; owned by connected bus

## Wiring

Port  | Type | Description
------|------|------------
`bus` | Bus  | Connected bus that owns terminal voltage variables and power-balance constraints

### Configuration Validation

Each load `id` must be nonempty and unique across all devices. Its `bus` value
must identify an existing bus.

## Model Constraints

### Internal Constraints

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

### External Constraints

Symbol           | Units  | Description                          | Note
-----------------|--------|--------------------------------------|------
$\mathrm{DIVP}$ | [p.u.] | Connected-bus active-power balance   | Enum name: `DIVP`; owned by connected bus
$\mathrm{DIVQ}$ | [p.u.] | Connected-bus reactive-power balance | Enum name: `DIVQ`; owned by connected bus

The fixed signed injection is added to the connected-bus balances:

```math
\begin{aligned}
  \mathrm{DIVP} &\mathrel{+}= u p \\
  \mathrm{DIVQ} &\mathrel{+}= u q
\end{aligned}
```

## Initialization

The Load model has no internal variables to initialize. Its fixed `p`, `q`,
and `online` values are read from the companion state.

## Model Outputs

Output | Units  | Description                                | Note
-------|--------|--------------------------------------------|------
`p`    | [p.u.] | Fixed signed active-power injection        | Preserved in the load device state
`q`    | [p.u.] | Fixed signed reactive-power injection      | Preserved in the load device state
`ir`   | [p.u.] | Bus-current injection, real component      | Written under the load `id` at the connected bus
`ii`   | [p.u.] | Bus-current injection, imaginary component | Written under the load `id` at the connected bus

The current injected into the bus follows $S=VI^{\ast}$:

```math
\begin{aligned}
  I_r &= u\dfrac{p\cos V_A + q\sin V_A}{V_M} \\
  I_i &= u\dfrac{p\sin V_A - q\cos V_A}{V_M}
\end{aligned}
```

For an online load, state output requires $V_M>0$. An offline load writes zero
current without evaluating the division above. Existing status fields and
unrelated state data are preserved.
