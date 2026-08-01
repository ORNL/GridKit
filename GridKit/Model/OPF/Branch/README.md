# **Branch Model**

The Branch model represents a two-terminal phasor-domain $\pi$ branch with an
optional off-nominal tap magnitude and phase shift on the from-bus side.
Terminal power contributions are oriented entering the adjacent buses. An
optional apparent-power rating adds one squared-flow constraint at each
terminal.

Notes:
- Setting $\tau = 1$ and $\theta = 0$ gives the ordinary symmetric transmission-line
  $\pi$ model.
- $G$ and $B$ are total branch shunt values split equally between terminals.
- The branch has no solver-owned variables; it contributes directly to the
  connected-bus power balances.
- An open branch remains in the model with zero-valued power and derivative
  contributions.

## Circuit Diagram

An ideal complex tap is placed on the from-bus side of the branch equivalent.
The ordinary transmission-line $\pi$ model is recovered with $\tau = 1$ and
$\theta = 0$.

![](../../../../docs/Figures/transformer-branch.png)

Figure 1: Branch equivalent circuit

## Model Parameters

Symbol        | Units  | Description                              | Typical Value | Note
--------------|--------|------------------------------------------|---------------|------
$R$           | [p.u.] | Branch series resistance                 |               | Input name: `R`
$X$           | [p.u.] | Branch series reactance                  |               | Input name: `X`
$G$           | [p.u.] | Total branch shunt conductance           | 0             | Input name: `G`
$B$           | [p.u.] | Total branch shunt susceptance           | 0             | Input name: `B`
$S^{\max}$   | [p.u.] | Apparent-power rating at each terminal   |               | Optional; input name: `smax`

## Model Inputs

Input   | Symbol | Units  | Description                              | Default
--------|--------|--------|------------------------------------------|--------
`tap`   | $\tau$ | [p.u.] | Off-nominal tap magnitude on from-bus side | 1
`phase` | $\theta$ | [rad] | Phase-shift angle                        | 0
`open`  | $o$    | [-]    | Branch open status                       | 0

These inputs are constant while a solve is running.

### Configuration Validation

Invalid Branch parameter and operating data are rejected by the following
checks:

```math
\begin{aligned}
  &R, X, G, B \in \mathbb{R}\ \text{and finite} \\
  &R^2 + X^2 \in \mathbb{R}\ \text{and finite},
    \qquad R^2 + X^2 > 0 \\
  &S^{\max} \in \mathbb{R}\ \text{and finite}
    \quad\text{when } S^{\max} \text{ is supplied} \\
  &S^{\max} \ge 0
    \quad\text{when } S^{\max} \text{ is supplied} \\
  &\left(S^{\max}\right)^2 \in \mathbb{R}\ \text{and finite}
    \quad\text{when } S^{\max} \text{ is supplied} \\
  &\tau > 0 \quad\text{and finite} \\
  &\theta \in \mathbb{R}\ \text{and finite}
\end{aligned}
```

Each branch `id` must be nonempty and unique across all devices. Its `from` and
`to` values must identify two distinct existing buses. The closed-branch
topology must connect all buses.

### Model Derived Parameters

The series and shunt admittances are:

```math
\begin{aligned}
  Y_{\mathrm{br}} &= \dfrac{1}{R + jX} \\
  Y_{\mathrm{sh}} &= G + jB
\end{aligned}
```

The terminal currents and powers are oriented entering the adjacent buses. The
nominal $\pi$-branch admittance matrix is:

```math
\begin{aligned}
  \mathbf{Y}_0
    &=
    \begin{bmatrix}
      -Y_{\mathrm{br}} & Y_{\mathrm{br}} \\
      Y_{\mathrm{br}}
      & -Y_{\mathrm{br}}
    \end{bmatrix}
    +
    \dfrac{1}{2}
    \begin{bmatrix}
      -Y_{\mathrm{sh}} & 0 \\
      0 & -Y_{\mathrm{sh}}
    \end{bmatrix}
\end{aligned}
```

The off-nominal transformer transformation uses the from bus as the tap side:

```math
\begin{aligned}
  \mathbf{M}
    &=
    \begin{bmatrix}
      \tau^{-1} & 0 \\
      0 & e^{j\theta}
    \end{bmatrix} \\
  \mathbf{Y} &= \mathbf{M}^{\dagger}\mathbf{Y}_0\mathbf{M}
\end{aligned}
```

Therefore:

```math
\begin{aligned}
  Y_{ff} &= -\dfrac{Y_{\mathrm{br}}+Y_{\mathrm{sh}}/2}{\tau^2} \\
  Y_{ft} &=  \dfrac{Y_{\mathrm{br}}e^{j\theta}}{\tau} \\
  Y_{tf} &=  \dfrac{Y_{\mathrm{br}}e^{-j\theta}}{\tau} \\
  Y_{tt} &= -\left(Y_{\mathrm{br}}+Y_{\mathrm{sh}}/2\right)
\end{aligned}
```

All four derived terminal-admittance entries must be finite.

For the constraints below, write $Y_{mn}=G_{mn}+jB_{mn}$ and define the
closed indicator as $c(o)=1$ when `open` is false or omitted and $c(o)=0$
when `open` is true.

## Model Variables

### Internal Variables

Symbol | Units | Description | Note
-------|-------|-------------|-----
None.  |       |             |

### External Variables

Symbol  | Units  | Description                       | Note
--------|--------|-----------------------------------|------
$V_{Mf}$ | [p.u.] | From-bus voltage magnitude       | Enum name: `VMF`; owned by from bus
$V_{Af}$ | [rad]  | From-bus voltage angle           | Enum name: `VAF`; owned by from bus
$V_{Mt}$ | [p.u.] | To-bus voltage magnitude         | Enum name: `VMT`; owned by to bus
$V_{At}$ | [rad]  | To-bus voltage angle             | Enum name: `VAT`; owned by to bus

## Wiring

Port   | Type | Description
-------|------|------------
`from` | Bus  | From bus that owns terminal voltage variables and power-balance constraints
`to`   | Bus  | To bus that owns terminal voltage variables and power-balance constraints

## Model Constraints

For readability, define $\delta=V_{Af}-V_{At}$ and:

```math
\begin{aligned}
  A &= G_{ft}\cos\delta + B_{ft}\sin\delta \\
  C &= G_{ft}\sin\delta - B_{ft}\cos\delta \\
  D &= G_{tf}\cos\delta - B_{tf}\sin\delta \\
  E &= -G_{tf}\sin\delta - B_{tf}\cos\delta
\end{aligned}
```

The terminal active- and reactive-power injections are:

```math
\begin{aligned}
  P_f &= c(o)\left(G_{ff}V_{Mf}^2 + V_{Mf}V_{Mt}A\right) \\
  Q_f &= c(o)\left(-B_{ff}V_{Mf}^2 + V_{Mf}V_{Mt}C\right) \\
  P_t &= c(o)\left(G_{tt}V_{Mt}^2 + V_{Mf}V_{Mt}D\right) \\
  Q_t &= c(o)\left(-B_{tt}V_{Mt}^2 + V_{Mf}V_{Mt}E\right)
\end{aligned}
```

### Internal Constraints

Symbol  | Units       | Description                              | Note
--------|-------------|------------------------------------------|------
$S_f^2$ | [p.u.^2] | Squared apparent power at from terminal  | Enum name: `SF2`; present only when `smax` is supplied
$S_t^2$ | [p.u.^2] | Squared apparent power at to terminal    | Enum name: `ST2`; present only when `smax` is supplied

When `smax` is supplied, the Branch model owns two bounded squared-flow
constraints:

```math
\begin{aligned}
  S_f^2 &= P_f^2 + Q_f^2,
  &0 \le S_f^2 \le \left(S^{\max}\right)^2 \\
  S_t^2 &= P_t^2 + Q_t^2,
  &0 \le S_t^2 \le \left(S^{\max}\right)^2
\end{aligned}
```

When `smax` is omitted, the Branch model owns no internal constraints.

### External Constraints

Symbol              | Units  | Description                         | Note
--------------------|--------|-------------------------------------|------
$\mathrm{DIVP}_f$  | [p.u.] | From-bus active-power balance       | Enum name: `DIVPF`; owned by from bus
$\mathrm{DIVQ}_f$  | [p.u.] | From-bus reactive-power balance     | Enum name: `DIVQF`; owned by from bus
$\mathrm{DIVP}_t$  | [p.u.] | To-bus active-power balance         | Enum name: `DIVPT`; owned by to bus
$\mathrm{DIVQ}_t$  | [p.u.] | To-bus reactive-power balance       | Enum name: `DIVQT`; owned by to bus

The terminal powers are added to the adjacent bus balances with positive sign
because they are oriented entering the buses:

```math
\begin{aligned}
  \mathrm{DIVP}_f &\mathrel{+}= P_f \\
  \mathrm{DIVQ}_f &\mathrel{+}= Q_f \\
  \mathrm{DIVP}_t &\mathrel{+}= P_t \\
  \mathrm{DIVQ}_t &\mathrel{+}= Q_t
\end{aligned}
```

## Initialization

The Branch model has no internal variables to initialize. The fixed `tap`,
`phase`, and `open` operating inputs use the defaults in the Model Inputs table
when their companion-state values are omitted.

## Model Outputs

Output  | Units  | Description                    | Note
--------|--------|--------------------------------|------
`tap`   | [p.u.] | Off-nominal tap magnitude      | Preserved from the input state
`phase` | [rad]  | Phase-shift angle              | Preserved from the input state
`open`  | [-]    | Branch open status             | Preserved from the input state

The Branch model does not write solved terminal-power fields. The copied
solution state preserves its operating inputs and unrelated state data.
