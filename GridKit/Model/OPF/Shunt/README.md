# **Shunt Model**

The Shunt model represents a fixed admittance connected from one bus to ground.
It owns no optimization variables or constraints and adds its voltage-dependent
active- and reactive-power injection directly to the connected-bus balances.

Notes:
- $G$ and $B$ use the network injection convention. Positive capacitive $B$
  injects positive reactive power into the connected bus.
- An offline shunt contributes zero power and current while retaining its
  admittance parameters.

## Model Parameters

Symbol | Units  | Description             | Typical Value | Note
-------|--------|-------------------------|---------------|------
$G$    | [p.u.] | Shunt conductance       | 0             | Input name: `G`
$B$    | [p.u.] | Shunt susceptance       | 0             | Input name: `B`

### Parameter Validation

Invalid Shunt parameter data are rejected by the following checks:

```math
G, B \in \mathbb{R}\ \text{and finite}
```

## Model Inputs

Input    | Symbol | Units | Description                 | Default
---------|--------|-------|-----------------------------|--------
`online` | $u$    | [-]   | Shunt connection status     | 1

The `online` input is constant while a solve is running.

### Model Derived Parameters

Define the shunt admittance and connection indicator as:

```math
\begin{aligned}
  Y_{\mathrm{sh}} &= G + jB \\
  u &=
  \begin{cases}
    1, & \text{when online is true or omitted} \\
    0, & \text{when online is false}
  \end{cases}
\end{aligned}
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

Each shunt `id` must be nonempty and unique across all devices. Its `bus` value
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

The shunt current is oriented entering the connected bus and is
$I_{\mathrm{sh}}=-uY_{\mathrm{sh}}V$. The corresponding power injection is:

```math
\begin{aligned}
  P_{\mathrm{sh}} &= -uG V_M^2 \\
  Q_{\mathrm{sh}} &=  uB V_M^2
\end{aligned}
```

The shunt adds these values to the connected-bus balances:

```math
\begin{aligned}
  \mathrm{DIVP} &\mathrel{+}= -uG V_M^2 \\
  \mathrm{DIVQ} &\mathrel{+}=  uB V_M^2
\end{aligned}
```

## Initialization

The Shunt model has no internal variables to initialize. Its contribution is
evaluated from the connected-bus voltage and the fixed `online` input.

## Model Outputs

Output | Units  | Description                                | Note
-------|--------|--------------------------------------------|------
`ir`   | [p.u.] | Bus-current injection, real component      | Written under the shunt `id` at the connected bus
`ii`   | [p.u.] | Bus-current injection, imaginary component | Written under the shunt `id` at the connected bus

For the accepted connected-bus voltage, define:

```math
\begin{aligned}
  V_r &= V_M \cos V_A \\
  V_i &= V_M \sin V_A
\end{aligned}
```

The shunt current injected into the bus is:

```math
\begin{aligned}
  I_r &= -u\left(GV_r - BV_i\right) \\
  I_i &= -u\left(BV_r + GV_i\right)
\end{aligned}
```

The model does not overwrite shunt device `p` or `q` fields. Existing status
fields and unrelated state data are preserved.
