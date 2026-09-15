# Machine Model

`Machine` represents a three-phase round-rotor synchronous machine with one
field winding and three damper windings. Internal quantities use machine
per-unit bases; terminal voltages and injected currents use SI units.

## Block Diagram

![Machine model block diagram](../../../../../../docs/Figures/EMT/Machine/diagram.png)

Figure 1: Machine model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$S_\mathrm{b}$ | [VA] | `S` | Rated three-phase apparent power | Required, positive
$V_\mathrm{b}$ | [V] | `V` | Rated line-to-line RMS voltage | Required, positive
$f_\mathrm{b}$ | [Hz] | `f` | Rated electrical frequency | Required, positive
$H$ | [s] | `H` | Inertia constant | Required, positive
$F$ | [p.u.] | `F` | Friction torque factor | Nonnegative
$R_s$ | [p.u.] | `Rs` | Stator winding resistance | Nonnegative
$L_l$ | [p.u.] | `Ll` | Stator leakage inductance | Positive
$L_\mathrm{md}$ | [p.u.] | `Lmd` | Unsaturated d-axis magnetizing inductance | Positive
$L_\mathrm{mq}$ | [p.u.] | `Lmq` | Unsaturated q-axis magnetizing inductance | Positive
$L_0$ | [p.u.] | `L0` | Zero-sequence inductance | Positive, defaults to $L_l$
$R_\mathrm{fd}$ | [p.u.] | `Rfd` | Field winding resistance | Positive
$L_\mathrm{lfd}$ | [p.u.] | `Llfd` | Field leakage inductance | Positive
$R_\mathrm{1d}$ | [p.u.] | `R1d` | d-axis damper resistance | Positive
$L_\mathrm{l1d}$ | [p.u.] | `Ll1d` | d-axis damper leakage inductance | Positive
$R_\mathrm{1q}$ | [p.u.] | `R1q` | q-axis damper 1 resistance | Positive
$L_\mathrm{l1q}$ | [p.u.] | `Ll1q` | q-axis damper 1 leakage inductance | Positive
$R_\mathrm{2q}$ | [p.u.] | `R2q` | q-axis damper 2 resistance | Positive
$L_\mathrm{l2q}$ | [p.u.] | `Ll2q` | q-axis damper 2 leakage inductance | Positive
$S(1.0)$ | [p.u.] | `S10` | Saturation factor at $1.0$ per-unit flux | Nonnegative
$S(1.2)$ | [p.u.] | `S12` | Saturation factor at $1.2$ per-unit flux | Nonnegative

### Parameter Validation

All parameters and derived coefficients must be finite.

```math
\begin{aligned}
S_\mathrm{b} &> 0 \\
V_\mathrm{b} &> 0 \\
f_\mathrm{b} &> 0 \\
H &> 0 \\
F &\ge 0 \\
R_s &\ge 0 \\
L_l, L_\mathrm{md}, L_\mathrm{mq}, L_0 &> 0 \\
L_\mathrm{lfd}, L_\mathrm{l1d}, L_\mathrm{l1q}, L_\mathrm{l2q} &> 0 \\
R_\mathrm{fd}, R_\mathrm{1d}, R_\mathrm{1q}, R_\mathrm{2q} &> 0 \\
S(1.0) = S(1.2) &= 0 \quad\text{or}\quad 0 \le S(1.0) < S(1.2)
\end{aligned}
```

### Derived Parameters

The phase count is $N=3$, with $\mathcal{N}=\{1,2,3\}$ in $(a,b,c)$ order.

```math
\begin{aligned}
\omega_\mathrm{b} &= 2\pi f_\mathrm{b} \\
V_\mathrm{pk} &= \sqrt{2/3}\,V_\mathrm{b} \\
I_\mathrm{pk} &= \dfrac{\sqrt{2}\,S_\mathrm{b}}{\sqrt{3}\,V_\mathrm{b}} \\
k_\mathrm{fd} &= \dfrac{R_\mathrm{fd}}{L_\mathrm{md}}
\end{aligned}
```

For $S(1.2)>0$,

```math
\begin{aligned}
s &= \sqrt{\dfrac{S(1.0)}{S(1.2)}} \\
S_A &= \dfrac{1.2s-1}{s-1} \\
S_B &= \dfrac{S(1.2)}{(S_A-1.2)^2}
\end{aligned}
```

Otherwise, $S_A=S_B=0$.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Bus voltage at machine port | $\mathbf{v} \in \mathbb{R}^N$
$\mathbf{i}$ | `i` | Output | [A] | Current injection at machine port | $\mathbf{i} \in \mathbb{R}^N$
$E_\mathrm{fd}$ | `efd` | Input | [p.u.] | Field voltage on the exciter base | Held constant when unattached
$P_m$ | `pm` | Input | [p.u.] | Mechanical power | Held constant when unattached
$\omega_r$ | `speed` | Output | [p.u.] | Rotor speed |

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\theta$ | [rad] | Electrical rotor angle | d-axis relative to the phase-a magnetic axis
$\omega_r$ | [p.u.] | Rotor speed | $\omega_r>0$
$\psi_d$ | [p.u.] | d-axis stator flux linkage |
$\psi_q$ | [p.u.] | q-axis stator flux linkage |
$\psi_0$ | [p.u.] | Zero-sequence stator flux linkage |
$\psi_\mathrm{fd}$ | [p.u.] | Field flux linkage |
$\psi_\mathrm{1d}$ | [p.u.] | d-axis damper flux linkage |
$\psi_\mathrm{1q}$ | [p.u.] | q-axis damper 1 flux linkage |
$\psi_\mathrm{2q}$ | [p.u.] | q-axis damper 2 flux linkage |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$i_d$ | [p.u.] | d-axis stator current | Positive out of the machine
$i_q$ | [p.u.] | q-axis stator current | Positive out of the machine
$i_0$ | [p.u.] | Zero-sequence stator current |
$i_\mathrm{fd}$ | [p.u.] | Field current |
$i_\mathrm{1d}$ | [p.u.] | d-axis damper current |
$i_\mathrm{1q}$ | [p.u.] | q-axis damper 1 current |
$i_\mathrm{2q}$ | [p.u.] | q-axis damper 2 current |
$\psi_\mathrm{ad}$ | [p.u.] | d-axis air-gap flux linkage |
$\psi_\mathrm{aq}$ | [p.u.] | q-axis air-gap flux linkage |
$\psi_\mathrm{at}$ | [p.u.] | Air-gap flux magnitude | $\psi_\mathrm{at}>0$
$K_s$ | [-] | Saturation factor |
$T_e$ | [p.u.] | Electrical torque |
$\mathbf{i}_s$ | [p.u.] | Instantaneous stator currents | $\mathbf{i}_s \in \mathbb{R}^N$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

Bus voltage is algebraic when its current-balance row has no voltage derivative.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$E_\mathrm{fd}$ | [p.u.] | Field voltage input |
$P_m$ | [p.u.] | Mechanical power input |

## Model Equations

The amplitude-invariant transformation uses
$\boldsymbol{\gamma}=[0,-2\pi/3,2\pi/3]^\mathsf{T}$:

```math
\mathbf{T}(\theta)=\dfrac{2}{3}
\begin{bmatrix}
\cos(\theta+\gamma_1) & \cos(\theta+\gamma_2) & \cos(\theta+\gamma_3) \\
-\sin(\theta+\gamma_1) & -\sin(\theta+\gamma_2) & -\sin(\theta+\gamma_3) \\
\frac{1}{2} & \frac{1}{2} & \frac{1}{2}
\end{bmatrix}
```

```math
\begin{aligned}
\begin{bmatrix}v_d & v_q & v_0\end{bmatrix}^\mathsf{T}
  &= \dfrac{1}{V_\mathrm{pk}}\mathbf{T}(\theta)\mathbf{v} \\
e_\mathrm{fd} &= k_\mathrm{fd}E_\mathrm{fd} \\
L_\mathrm{ad} &= K_sL_\mathrm{md} \\
L_\mathrm{aq} &= K_sL_\mathrm{mq}
\end{aligned}
```

The common saturation factor uses the
[quadratic ramp](../../../../../CommonMath.md#quadratic-ramp) $q$.
Differential leakage between rotor windings is neglected.

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= \dfrac{\mathrm{d}\theta}{\mathrm{d}t} - \omega_\mathrm{b}\,\omega_r \\
0 &= 2H\dfrac{\mathrm{d}\omega_r}{\mathrm{d}t}
     - \dfrac{P_m}{\omega_r} + T_e + F\,\omega_r \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_d}{\mathrm{d}t}
     - \omega_r\,\psi_q - R_s\,i_d - v_d \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_q}{\mathrm{d}t}
     + \omega_r\,\psi_d - R_s\,i_q - v_q \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_0}{\mathrm{d}t}
     - R_s\,i_0 - v_0 \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_\mathrm{fd}}{\mathrm{d}t}
     + R_\mathrm{fd}\,i_\mathrm{fd} - e_\mathrm{fd} \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_\mathrm{1d}}{\mathrm{d}t}
     + R_\mathrm{1d}\,i_\mathrm{1d} \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_\mathrm{1q}}{\mathrm{d}t}
     + R_\mathrm{1q}\,i_\mathrm{1q} \\
0 &= \dfrac{1}{\omega_\mathrm{b}}\dfrac{\mathrm{d}\psi_\mathrm{2q}}{\mathrm{d}t}
     + R_\mathrm{2q}\,i_\mathrm{2q}
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
0 &= \psi_d + (L_\mathrm{ad} + L_l)\,i_d
     - L_\mathrm{ad}\,(i_\mathrm{fd} + i_\mathrm{1d}) \\
0 &= \psi_q + (L_\mathrm{aq} + L_l)\,i_q
     - L_\mathrm{aq}\,(i_\mathrm{1q} + i_\mathrm{2q}) \\
0 &= \psi_0 + L_0\,i_0 \\
0 &= \psi_\mathrm{fd} - (L_\mathrm{ad} + L_\mathrm{lfd})\,i_\mathrm{fd}
     - L_\mathrm{ad}\,i_\mathrm{1d} + L_\mathrm{ad}\,i_d \\
0 &= \psi_\mathrm{1d} - L_\mathrm{ad}\,i_\mathrm{fd}
     - (L_\mathrm{ad} + L_\mathrm{l1d})\,i_\mathrm{1d} + L_\mathrm{ad}\,i_d \\
0 &= \psi_\mathrm{1q} - (L_\mathrm{aq} + L_\mathrm{l1q})\,i_\mathrm{1q}
     - L_\mathrm{aq}\,i_\mathrm{2q} + L_\mathrm{aq}\,i_q \\
0 &= \psi_\mathrm{2q} - L_\mathrm{aq}\,i_\mathrm{1q}
     - (L_\mathrm{aq} + L_\mathrm{l2q})\,i_\mathrm{2q} + L_\mathrm{aq}\,i_q \\
0 &= \psi_\mathrm{ad} - \psi_d - L_l\,i_d \\
0 &= \psi_\mathrm{aq} - \psi_q - L_l\,i_q \\
0 &= \psi_\mathrm{at} - \sqrt{\psi_\mathrm{ad}^2 + \psi_\mathrm{aq}^2} \\
0 &= K_s\left(1 + S_B\,q(\psi_\mathrm{at} - S_A)\right) - 1 \\
0 &= T_e - (\psi_d\,i_q - \psi_q\,i_d) \\
0 &= i_{s,n} - \left(i_d \cos(\theta + \gamma_n)
     - i_q \sin(\theta + \gamma_n) + i_0\right),
     \quad n \in \mathcal{N}
\end{aligned}
```

### External Equations

```math
\mathbf{i} \leftarrow I_\mathrm{pk}\,\mathbf{i}_s
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`theta` | [rad] | Electrical rotor angle |
`omega` | [p.u.] | Rotor speed |
`te` | [p.u.] | Electrical torque |
`ifd` | [p.u.] | Field current |
`efd` | [p.u.] | Applied field voltage | Exciter base
`ks` | [-] | Saturation factor |
`psi_at` | [p.u.] | Air-gap flux magnitude |
`ia` | [A] | Phase a current injection |
`ib` | [A] | Phase b current injection |
`ic` | [A] | Phase c current injection |
`p` | [W] | Active power injection | $\mathbf{v}^\mathsf{T}\mathbf{i}$
`q` | [var] | Reactive power injection | $S_\mathrm{b}(v_qi_d-v_di_q)$
`id` | [p.u.] | d-axis stator current |
`iq` | [p.u.] | q-axis stator current |
