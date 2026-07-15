# **Renewable Energy Electrical Control Model (REECA)**

REECA is a WECC renewable electrical-control model for inverter-coupled
resources.

## Notes

- When used with REPCA, connect REPCA `qext` to REECA `qext` and REPCA `pext`
  to REECA `pref`.
- The `pe` and `qgen` feedback signal ports must either both be connected or
  both be omitted. When omitted, initialization derives constant feedback from
  the initialized `ipcmd` and `iqcmd` signal starts.
- The `omegag` signal input is required when `PFlag` is 1. When `PFlag` is 0,
  an omitted `omegag` input defaults to zero. Optional `qext`, `pfaref`, and
  `pref` inputs default to their initialized constant values when omitted.
- Timer-based post-dip reactive-current injection hold and active-current limit
  hold are not modeled in this version. `Thld` and `Thld2` must be zero, and
  `Iqfrz` is therefore unused.

## Block Diagram

Standard REECA block diagram.

![](../../../../../docs/Figures/PhasorDynamics/REECA/diagram.png)

Figure 1: REECA block diagram. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol                              | Units    | JSON     | Description                                             | Typical Value | Note
------------------------------------|----------|----------|---------------------------------------------------------|---------------|------
$S^\mathrm{base}$                   | [MVA]    | `mva`    | REECA component power base                              | 100.0         | Required positive value; block name: `MVABase`
$s_{\mathrm{pf}}$                   | [binary] | `PfFlag` | Power-factor control flag                               | 0             | Block name: `PfFlag`; 1 = power-factor control, 0 = Q control
$s_V$                               | [binary] | `VFlag`  | Voltage-control mode flag                               | 0             | Block name: `VFlag`; 1 = Q control, 0 = voltage control
$s_Q$                               | [binary] | `QFlag`  | Reactive-power control flag                             | 0             | Block name: `QFlag`; 1 = voltage/Q control, 0 = constant pf or Q control
$s_P$                               | [binary] | `PFlag`  | Active-power reference speed-proxy multiplier flag      | 0             | Block name: `PFlag`; 0 = direct reference, 1 = multiply by internal speed proxy
$s_{PQ}$                            | [binary] | `Pqflag` | P/Q priority flag for converter current limit           | 0             | Block name: `Pqflag`; 0 = Q priority, 1 = P priority
$T_{\mathrm{rv}}$                   | [sec]    | `Trv`    | Voltage-measurement filter time constant                | 0.02          | State 1
$T_{\mathrm{p}}$                    | [sec]    | `Tp`     | Electrical-power measurement filter time constant       | 0.0           | State 2
$V_0^\mathrm{ref}$                  | [p.u.]   | `Vref0`  | Outer-loop voltage reference                            | $V_{T,0}$     | Initialized from terminal voltage if omitted
$V_{\mathrm{dip}}$                  | [p.u.]   | `Vdip`   | Low-voltage threshold for the voltage-band gate         | 0.85          |
$V_{\mathrm{up}}$                   | [p.u.]   | `Vup`    | High-voltage threshold for the voltage-band gate        | 1.15          |
$D_1^\mathrm{db}$                   | [p.u.]   | `dbd1`   | Lower deadband threshold for voltage-error response     | 0.0           |
$D_2^\mathrm{db}$                   | [p.u.]   | `dbd2`   | Upper deadband threshold for voltage-error response     | 0.0           |
$K_{\mathrm{qv}}$                   | [p.u.]   | `kqv`    | Reactive-current injection gain outside the voltage band | 5.0          |
$I_{q,\mathrm{inj}}^{\min}$         | [p.u.]   | `Iql1`   | Minimum reactive-current injection limit                | -1.1          |
$I_{q,\mathrm{inj}}^{\max}$         | [p.u.]   | `Iqh1`   | Maximum reactive-current injection limit                | 1.1           |
$I_{q,\mathrm{inj}}^\mathrm{frz}$   | [p.u.]   | `Iqfrz`  | Held reactive-current injection value after voltage dip | 0.0           | Unused when $T_{\mathrm{hld}} = 0$
$T_{\mathrm{hld}}$                  | [sec]    | `Thld`   | Reactive-current injection hold time after voltage dip clears | 0.0      | Required to be zero in this version
$T_{\mathrm{hld2}}$                 | [sec]    | `Thld2`  | Active-current limit hold time after voltage dip clears | 0.0           | Optional; required to be zero in this version
$Q^{\max}$                          | [p.u.]   | `Qmax`   | Maximum reactive-power control limit                    | 0.436         |
$Q^{\min}$                          | [p.u.]   | `Qmin`   | Minimum reactive-power control limit                    | -0.436        |
$K_{\mathrm{qp}}$                   | [p.u.]   | `Kqp`    | Reactive-power control proportional gain                | 0.0           |
$K_{\mathrm{qi}}$                   | [p.u./s] | `Kqi`    | Reactive-power control integral gain                    | 0.1           |
$V^{\max}$                          | [p.u.]   | `Vmax`   | Maximum voltage-control limit                           | 1.1           |
$V^{\min}$                          | [p.u.]   | `Vmin`   | Minimum voltage-control limit                           | 0.9           |
$V_1^\mathrm{ref}$                  | [p.u.]   | `Vref1`  | Inner-loop voltage-control reference/bias               | 0.0           |
$K_{\mathrm{vp}}$                   | [p.u.]   | `Kvp`    | Voltage-control proportional gain                       | 18.0          |
$K_{\mathrm{vi}}$                   | [p.u./s] | `Kvi`    | Voltage-control integral gain                           | 5.0           |
$T_{\mathrm{iq}}$                   | [sec]    | `Tiq`    | Reactive-current command lag time constant              | 0.02          | State 5
$T_{\mathrm{pord}}$                 | [sec]    | `Tpord`  | Active-power order filter time constant                 | 0.02          | State 6
$R_P^{\max}$                        | [p.u./s] | `dPmax`  | Positive active-power order ramp-rate limit             | 99.0          |
$R_P^{\min}$                        | [p.u./s] | `dPmin`  | Negative active-power order ramp-rate limit             | -99.0         |
$P^{\max}$                          | [p.u.]   | `Pmax`   | Maximum active-power order limit                        | 1.0           |
$P^{\min}$                          | [p.u.]   | `Pmin`   | Minimum active-power order limit                        | 0.0           |
$I^{\max}$                          | [p.u.]   | `Imax`   | Maximum total converter current                         | 1.3           |
$V_{q,1}$                           | [p.u.]   | `Vq1`    | VDL1 voltage point 1                                    | 0.2           | VDL1 in Fig. 1
$I_{q,1}^{\max}$                    | [p.u.]   | `Iq1`    | VDL1 reactive-current limit point 1                     | 1.3           | VDL1 in Fig. 1
$V_{q,2}$                           | [p.u.]   | `Vq2`    | VDL1 voltage point 2                                    | 0.5           | VDL1 in Fig. 1
$I_{q,2}^{\max}$                    | [p.u.]   | `Iq2`    | VDL1 reactive-current limit point 2                     | 1.3           | VDL1 in Fig. 1
$V_{q,3}$                           | [p.u.]   | `Vq3`    | VDL1 voltage point 3                                    | 0.75          | VDL1 in Fig. 1
$I_{q,3}^{\max}$                    | [p.u.]   | `Iq3`    | VDL1 reactive-current limit point 3                     | 1.3           | VDL1 in Fig. 1
$V_{q,4}$                           | [p.u.]   | `Vq4`    | VDL1 voltage point 4                                    | 1.0           | VDL1 in Fig. 1
$I_{q,4}^{\max}$                    | [p.u.]   | `Iq4`    | VDL1 reactive-current limit point 4                     | 1.3           | VDL1 in Fig. 1
$V_{p,1}$                           | [p.u.]   | `Vp1`    | VDL2 voltage point 1                                    | 0.2           | VDL2 in Fig. 1
$I_{p,1}^{\max}$                    | [p.u.]   | `Ip1`    | VDL2 active-current limit point 1                       | 1.3           | VDL2 in Fig. 1
$V_{p,2}$                           | [p.u.]   | `Vp2`    | VDL2 voltage point 2                                    | 0.5           | VDL2 in Fig. 1
$I_{p,2}^{\max}$                    | [p.u.]   | `Ip2`    | VDL2 active-current limit point 2                       | 1.3           | VDL2 in Fig. 1
$V_{p,3}$                           | [p.u.]   | `Vp3`    | VDL2 voltage point 3                                    | 0.75          | VDL2 in Fig. 1
$I_{p,3}^{\max}$                    | [p.u.]   | `Ip3`    | VDL2 active-current limit point 3                       | 1.3           | VDL2 in Fig. 1
$V_{p,4}$                           | [p.u.]   | `Vp4`    | VDL2 voltage point 4                                    | 1.0           | VDL2 in Fig. 1
$I_{p,4}^{\max}$                    | [p.u.]   | `Ip4`    | VDL2 active-current limit point 4                       | 1.3           | VDL2 in Fig. 1

### Parameter Validation

Invalid REECA parameter sets are rejected by the following checks. The displayed
equations use effective time constants with $\epsilon_T=10^{-3}$.

```math
\begin{aligned}
  T &\leftarrow \max\!\left(T, \epsilon_T\right)
    \quad T\in\{T_{\mathrm{rv}},T_{\mathrm{p}}\} \\
  S^\mathrm{base} &> 0 \\
  s_{\mathrm{pf}}, s_V, s_Q, s_P, s_{PQ}
    &\in \{0,1\} \\
  T_{\mathrm{rv}}, T_{\mathrm{p}}
    &\ge 0 \\
  0 \le V_{\mathrm{dip}}
    &< V_{\mathrm{up}} \\
  D_1^\mathrm{db}
    &\le 0 \le D_2^\mathrm{db} \\
  I_{q,\mathrm{inj}}^{\min}
    &\le I_{q,\mathrm{inj}}^{\max} \\
  T_{\mathrm{hld}}, T_{\mathrm{hld2}}
    &= 0 \\
  Q^{\min}
    &\le Q^{\max} \\
  V^{\min}
    &\le V^{\max} \\
  T_{\mathrm{iq}}, T_{\mathrm{pord}}
    &> 0 \\
  R_P^{\min}
    &< 0 < R_P^{\max} \\
  P^{\min}
    &\le P^{\max} \\
  I^{\max}
    &\ge 0 \\
  0 \le V_{q,1}
    &< V_{q,2} < V_{q,3} < V_{q,4} \\
  I_{q,k}^{\max}
    &\ge 0
    \quad k\in\{1,\ldots,4\} \\
  0 \le V_{p,1}
    &< V_{p,2} < V_{p,3} < V_{p,4} \\
  I_{p,k}^{\max}
    &\ge 0
    \quad k\in\{1,\ldots,4\}
\end{aligned}
```

### Model Derived Parameters

```math
\begin{aligned}
  s_{\mathrm{pf}}^\mathrm{off}
    &= 1 - s_{\mathrm{pf}} \\
  s_V^\mathrm{off}
    &= 1 - s_V \\
  s_Q^\mathrm{off}
    &= 1 - s_Q \\
  s_P^\mathrm{off}
    &= 1 - s_P \\
  s_{PQ}^\mathrm{off}
    &= 1 - s_{PQ} \\
  k_{\mathrm{base}}
    &= \dfrac{S^\mathrm{sys}}{S^\mathrm{base}}
\end{aligned}
```

```math
\begin{aligned}
  g_q(x)
    &=
      I_{q,1}^{\max}
      + \sum_{k\in\{1,\ldots,3\}}
        \text{linseg}\!\left(
          x;\,
          V_{q,k},\,
          V_{q,k+1},\,
          I_{q,k+1}^{\max} - I_{q,k}^{\max}
        \right) \\
  g_p(x)
    &=
      I_{p,1}^{\max}
      + \sum_{k\in\{1,\ldots,3\}}
        \text{linseg}\!\left(
          x;\,
          V_{p,k},\,
          V_{p,k+1},\,
          I_{p,k+1}^{\max} - I_{p,k}^{\max}
        \right) \\
  m_+(x,y)
    &= \dfrac{xy}{\text{max}(x,y)}
    \qquad x,y\ge 0
\end{aligned}
```

CommonMath defines the [linear segment](../../../../CommonMath.md#linseg)
helper used by $g_q$ and $g_p$. The VDL functions provide flat extrapolation
beyond the first and fourth voltage points. The smooth nonnegative current-limit
target $m_+$ uses the CommonMath [max](../../../../CommonMath.md#derived-functions)
helper. It is exactly zero when either operand is zero and is nonnegative for
the validated nonnegative VDL ordinates and initialized nonnegative
current-circle branches.

## Model Ports

Name     | Port   | Init    | Description
---------|--------|---------|------
`bus`    | Bus    | Known   | Terminal bus voltage
`pe`     | Input  | Unknown | Electrical active-power feedback
`qgen`   | Input  | Unknown | Reactive-power feedback
`omegag` | Input  | Known   | Internal speed proxy
`qext`   | Input  | Unknown | External reactive-power command
`pfaref` | Input  | Unknown | Power-factor angle reference
`pref`   | Input  | Unknown | External active-power reference
`iqcmd`  | Output | Known   | Reactive-current command output
`ipcmd`  | Output | Known   | Active-current command output

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                         | Note
------------------------|--------|-------------------------------------|------
$V^\mathrm{meas}$       | [p.u.] | Filtered terminal voltage           | State 1 in Fig. 1
$P^\mathrm{meas}$       | [p.u.] | Filtered electrical power           | State 2 in Fig. 1
$x_Q^\mathrm{PI}$       | [p.u.] | Reactive-power PI controller state  | State 3 in Fig. 1
$x_V^\mathrm{PI}$       | [p.u.] | Voltage PI controller state         | State 4 in Fig. 1
$Q_V$                   | [p.u.] | Reactive-current command lag state  | State 5 in Fig. 1
$P^\mathrm{ord}$        | [p.u.] | Filtered active-power order         | State 6 in Fig. 1

#### Algebraic

Symbol                              | Units    | Description                         | Note
------------------------------------|----------|-------------------------------------|------
$V_T$                               | [p.u.]   | Terminal voltage magnitude          |
$V_{\mathrm{safe}}^\mathrm{meas}$   | [p.u.]   | Safe filtered terminal voltage for divider blocks | Lower bounded by 0.01
$s_{\mathrm{dip}}$                  | [-]      | Voltage inside-band control gate    |
$e_V^\mathrm{db}$                   | [p.u.]   | Deadbanded voltage error            |
$I_q^\mathrm{inj}$                  | [p.u.]   | Reactive-current injection candidate | Component base
$Q^\mathrm{ref}$                    | [p.u.]   | Selected reactive-power reference   |
$e_Q$                               | [p.u.]   | Reactive-power control error        |
$V_Q^\mathrm{PI}$                   | [p.u.]   | Reactive-power control PI output    |
$e_V^\mathrm{PI}$                   | [p.u.]   | Voltage-control PI error            |
$f_P^\mathrm{ord}$                  | [p.u./s] | Active-power order derivative target | Before ramp-rate limit
$r_P^\mathrm{ord}$                  | [p.u./s] | Ramp-rate-limited active-power order derivative target |
$I_q^\mathrm{circ}$                 | [p.u.]   | Reactive-current limit from converter current circle |
$I_p^\mathrm{circ}$                 | [p.u.]   | Active-current limit from converter current circle |
$I_q^{\max}$                        | [p.u.]   | Final reactive-current upper limit  | Component base
$I_p^{\max}$                        | [p.u.]   | Final active-current upper limit    | Component base
$I_q^\mathrm{base}$                 | [p.u.]   | Base reactive-current command       | Component base
$I_q^\mathrm{raw}$                  | [p.u.]   | Raw reactive-current command before final limit | Component base
$I_q^\mathrm{cmd}$                  | [p.u.]   | Reactive-current command output     | System base
$I_p^\mathrm{cmd}$                  | [p.u.]   | Active-current command output       | System base

### External Variables

#### Differential
None.

#### Algebraic

Symbol                               | Units  | Type    | Description                       | Note
-------------------------------------|--------|---------|-----------------------------------|------
$V_{\mathrm{r}}$                     | [p.u.] | Known   | Terminal voltage, real component  | Bus input
$V_{\mathrm{i}}$                     | [p.u.] | Known   | Terminal voltage, imaginary component | Bus input
$P_e$                                | [p.u.] | Unknown | Electrical active-power feedback  | Signal port `pe`; system base
$Q^\mathrm{gen}$                     | [p.u.] | Unknown | Reactive-power feedback           | Signal port `qgen`; system base
$\omega_g$                           | [p.u.] | Known   | Internal speed proxy              | Signal port `omegag`; required when $s_P=1$, otherwise defaults to zero
$Q^\mathrm{ext}$                     | [p.u.] | Unknown | External reactive-power command   | Optional signal port `qext`; system base
$\phi^\mathrm{ref}$                  | [rad]  | Unknown | Power-factor angle reference      | Optional signal port `pfaref`
$P^\mathrm{ref}$                     | [p.u.] | Unknown | External active-power reference   | Optional signal port `pref`; system base

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{V}^\mathrm{meas}
    + \dfrac{1}{T_{\mathrm{rv}}}
      \left(V_T - V^\mathrm{meas}\right) \\
  0 &=
    -\dot{P}^\mathrm{meas}
    + \dfrac{1}{T_{\mathrm{p}}}
      \left(k_{\mathrm{base}}P_e - P^\mathrm{meas}\right) \\
  0 &=
    -\dot{x}_Q^\mathrm{PI}
    + s_{\mathrm{dip}}\,
      \text{antiwindup}
      \left(
        K_{\mathrm{qp}}e_Q + x_Q^\mathrm{PI},\,
        K_{\mathrm{qi}}e_Q;\,
        V^{\min}, V^{\max}
      \right) \\
  0 &=
    -\dot{x}_V^\mathrm{PI}
    + s_{\mathrm{dip}}\,
      \text{antiwindup}
      \left(
        K_{\mathrm{vp}}e_V^\mathrm{PI} + x_V^\mathrm{PI},\,
        K_{\mathrm{vi}}e_V^\mathrm{PI};\,
        -I_q^{\max}, I_q^{\max}
      \right) \\
  0 &=
    -\dot{Q}_V
    + \dfrac{s_{\mathrm{dip}}}{T_{\mathrm{iq}}}
      \left(
        \dfrac{Q^\mathrm{ref}}{V_{\mathrm{safe}}^\mathrm{meas}}
        - Q_V
      \right) \\
  0 &=
    -\dot{P}^\mathrm{ord}
    + s_{\mathrm{dip}}\,
      \text{antiwindup}
      \left(P^\mathrm{ord}, r_P^\mathrm{ord};\, P^{\min}, P^{\max}\right)
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#antiwindup)
target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &=
    -V_T^2
    + V_{\mathrm{r}}^2
    + V_{\mathrm{i}}^2 \\
  0 &=
    -V_{\mathrm{safe}}^\mathrm{meas}
    + \text{max}
      \left(V^\mathrm{meas}, 0.01\right) \\
  0 &=
    -s_{\mathrm{dip}}
    + \text{inside}
      \left(V_T;\, V_{\mathrm{dip}}, V_{\mathrm{up}}\right) \\
  0 &=
    -e_V^\mathrm{db}
    + \text{deadband2}
      \left(V_0^\mathrm{ref} - V^\mathrm{meas};\,
            D_1^\mathrm{db}, D_2^\mathrm{db}\right) \\
  0 &=
    -I_q^\mathrm{inj}
    + \text{clamp}
      \left(K_{\mathrm{qv}}e_V^\mathrm{db};\,
            I_{q,\mathrm{inj}}^{\min}, I_{q,\mathrm{inj}}^{\max}\right) \\
  0 &=
    -Q^\mathrm{ref}
    + s_{\mathrm{pf}}P^\mathrm{meas}\tan\!\left(\phi^\mathrm{ref}\right)
    + s_{\mathrm{pf}}^\mathrm{off}k_{\mathrm{base}}Q^\mathrm{ext} \\
  0 &=
    -e_Q
    + \text{clamp}
      \left(Q^\mathrm{ref};\, Q^{\min}, Q^{\max}\right)
    - k_{\mathrm{base}}Q^\mathrm{gen} \\
  0 &=
    -V_Q^\mathrm{PI}
    + \text{clamp}
      \left(K_{\mathrm{qp}}e_Q + x_Q^\mathrm{PI};\,
            V^{\min}, V^{\max}\right) \\
  0 &=
    -e_V^\mathrm{PI}
    + s_V V_Q^\mathrm{PI}
    + s_V^\mathrm{off}\left(Q^\mathrm{ref} + V_1^\mathrm{ref}\right)
    - V^\mathrm{meas} \\
  0 &=
    -f_P^\mathrm{ord}
    + \dfrac{1}{T_{\mathrm{pord}}}
      \left(\left(s_P^\mathrm{off} + s_P\omega_g\right)k_{\mathrm{base}}P^\mathrm{ref} - P^\mathrm{ord}\right) \\
  0 &=
    -r_P^\mathrm{ord}
    + \text{clamp}
      \left(f_P^\mathrm{ord};\, R_P^{\min}, R_P^{\max}\right) \\
  0 &=
    -\left(I_q^\mathrm{circ}\right)^2
    + \left(I^{\max}\right)^2
    - s_{PQ}\left(k_{\mathrm{base}}I_p^\mathrm{cmd}\right)^2 \\
  0 &=
    -\left(I_p^\mathrm{circ}\right)^2
    + \left(I^{\max}\right)^2
    - s_{PQ}^\mathrm{off}\left(k_{\mathrm{base}}I_q^\mathrm{cmd}\right)^2 \\
  0 &=
    -I_q^{\max}
    + m_+\left(g_q\!\left(V^\mathrm{meas}\right), I_q^\mathrm{circ}\right) \\
  0 &=
    -I_p^{\max}
    + m_+\left(g_p\!\left(V^\mathrm{meas}\right), I_p^\mathrm{circ}\right) \\
  0 &=
    -I_q^\mathrm{base}
    + \text{clamp}
      \left(K_{\mathrm{vp}}e_V^\mathrm{PI} + x_V^\mathrm{PI};\,
            -I_q^{\max}, I_q^{\max}\right) \\
  0 &=
    -I_q^\mathrm{raw}
    + s_Q I_q^\mathrm{base}
    + s_Q^\mathrm{off}Q_V
    + \left(1-s_{\mathrm{dip}}\right)I_q^\mathrm{inj} \\
  0 &=
    -k_{\mathrm{base}}I_q^\mathrm{cmd}
    + \text{clamp}
      \left(I_q^\mathrm{raw};\,
            -I_q^{\max}, I_q^{\max}\right) \\
  0 &=
    -k_{\mathrm{base}}I_p^\mathrm{cmd}
    + \text{clamp}
      \left(
        \dfrac{P^\mathrm{ord}}{V_{\mathrm{safe}}^\mathrm{meas}};\,
        0,\,
        I_p^{\max}
      \right)
\end{aligned}
```

Initialization selects the nonnegative branches of the squared $V_T$,
$I_q^\mathrm{circ}$, and $I_p^\mathrm{circ}$ algebraic residuals.

CommonMath defines helper targets and smooth approximations for
[max, clamp, deadband2, and inside](../../../../CommonMath.md#derived-functions).

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_{\mathrm{r}}, V_{\mathrm{i}}
    &\leftarrow \text{terminal-bus voltage} \\
  I_q^\mathrm{cmd}, I_p^\mathrm{cmd}
    &\leftarrow \text{current-command start} \\
  \omega_g
    &\leftarrow \text{internal-speed-proxy input, or }0
      \quad\text{when omitted and }s_P=0
\end{aligned}
```

### Internal Initialization

Define

```math
\begin{aligned}
  \text{awinit}(x^\star,f;\ell,u)
    &=
      \begin{cases}
        x^\star & f = 0 \\
        u + \epsilon_{\mathrm{sat}} & f > 0 \\
        \ell - \epsilon_{\mathrm{sat}} & f < 0
      \end{cases}
\end{aligned}
```

with $\epsilon_{\mathrm{sat}}>0$.

Initialization is performed by evaluating the steady-state residuals in
dependency order. Let subscript $0$ denote initial values and set all internal
derivatives to zero:

```math
\begin{aligned}
  V_{T,0}
    &= \sqrt{V_{\mathrm{r},0}^2 + V_{\mathrm{i},0}^2} \\
  V_0^\mathrm{ref}
    &= V_{T,0}\quad\text{if omitted} \\
  V_0^\mathrm{meas}
    &= V_{T,0} \\
  V_{\mathrm{safe},0}^\mathrm{meas}
    &= \text{max}\left(V_0^\mathrm{meas}, 0.01\right) \\
  P_0^\mathrm{meas}
    &= k_{\mathrm{base}}P_{e,0} \\
  s_{\mathrm{dip},0}
    &= \text{inside}
       \left(V_{T,0};\, V_{\mathrm{dip}}, V_{\mathrm{up}}\right) \\
  e_{V,0}^\mathrm{db}
    &=
      \text{deadband2}
      \left(V_0^\mathrm{ref} - V_0^\mathrm{meas};\,
            D_1^\mathrm{db}, D_2^\mathrm{db}\right) \\
  I_{q,0}^\mathrm{inj}
    &=
      \text{clamp}
      \left(K_{\mathrm{qv}}e_{V,0}^\mathrm{db};\,
            I_{q,\mathrm{inj}}^{\min}, I_{q,\mathrm{inj}}^{\max}\right) \\
  Q_0^\mathrm{ref}
    &=
      s_{\mathrm{pf}}P_0^\mathrm{meas}
      \tan\!\left(\phi_0^\mathrm{ref}\right)
      + s_{\mathrm{pf}}^\mathrm{off}k_{\mathrm{base}}Q_0^\mathrm{ext} \\
  e_{Q,0}
    &=
      \text{clamp}
      \left(Q_0^\mathrm{ref};\, Q^{\min}, Q^{\max}\right)
      - k_{\mathrm{base}}Q_0^\mathrm{gen} \\
  Q_{V,0}
    &= \dfrac{Q_0^\mathrm{ref}}{V_{\mathrm{safe},0}^\mathrm{meas}} \\
  P_0^\mathrm{ord}
    &= k_{\mathrm{base}}P_{e,0} \\
  f_{P,0}^\mathrm{ord}
    &= 0 \\
  r_{P,0}^\mathrm{ord}
    &= 0 \\
  u_{Q,0}^\mathrm{PI}
    &=
      \text{awinit}
      \left(
        s_V V_0^\mathrm{meas}
        + s_V^\mathrm{off}\left(Q_0^\mathrm{ref} + V_1^\mathrm{ref}\right),\,
        K_{\mathrm{qi}}e_{Q,0};\,
        V^{\min}, V^{\max}
      \right) \\
  V_{Q,0}^\mathrm{PI}
    &=
      \text{clamp}
      \left(u_{Q,0}^\mathrm{PI};\, V^{\min}, V^{\max}\right) \\
  e_{V,0}^\mathrm{PI}
    &=
      s_V V_{Q,0}^\mathrm{PI}
      + s_V^\mathrm{off}\left(Q_0^\mathrm{ref} + V_1^\mathrm{ref}\right)
      - V_0^\mathrm{meas} \\
  x_{Q,0}^\mathrm{PI}
    &= u_{Q,0}^\mathrm{PI} - K_{\mathrm{qp}}e_{Q,0}
\end{aligned}
```

The current limits and output commands are initialized from the final algebraic
equations:

```math
\begin{aligned}
  \left(I_{q,0}^\mathrm{circ}\right)^2
    &=
      \left(I^{\max}\right)^2
      - s_{PQ}\left(k_{\mathrm{base}}I_{p,0}^\mathrm{cmd}\right)^2 \\
  \left(I_{p,0}^\mathrm{circ}\right)^2
    &=
      \left(I^{\max}\right)^2
      - s_{PQ}^\mathrm{off}\left(k_{\mathrm{base}}I_{q,0}^\mathrm{cmd}\right)^2 \\
  I_{q,0}^{\max}
    &= m_+\left(g_q\!\left(V_0^\mathrm{meas}\right), I_{q,0}^\mathrm{circ}\right) \\
  I_{p,0}^{\max}
    &= m_+\left(g_p\!\left(V_0^\mathrm{meas}\right), I_{p,0}^\mathrm{circ}\right) \\
  u_{V,0}^\mathrm{PI}
    &=
      \text{awinit}
      \left(
        \dfrac{k_{\mathrm{base}}Q_0^\mathrm{gen}}
              {V_{\mathrm{safe},0}^\mathrm{meas}},\,
        K_{\mathrm{vi}}e_{V,0}^\mathrm{PI};\,
        -I_{q,0}^{\max}, I_{q,0}^{\max}
      \right) \\
  I_{q,0}^\mathrm{base}
    &=
      \text{clamp}
      \left(u_{V,0}^\mathrm{PI};\,
            -I_{q,0}^{\max},
            I_{q,0}^{\max}\right) \\
  x_{V,0}^\mathrm{PI}
    &= u_{V,0}^\mathrm{PI} - K_{\mathrm{vp}}e_{V,0}^\mathrm{PI} \\
  I_{q,0}^\mathrm{raw}
    &=
      s_Q I_{q,0}^\mathrm{base}
      + s_Q^\mathrm{off}Q_{V,0}
      + \left(1-s_{\mathrm{dip},0}\right)I_{q,0}^\mathrm{inj} \\
  k_{\mathrm{base}}I_{q,0}^\mathrm{cmd}
    &=
      \text{clamp}
      \left(I_{q,0}^\mathrm{raw};\,
            -I_{q,0}^{\max}, I_{q,0}^{\max}\right) \\
  k_{\mathrm{base}}I_{p,0}^\mathrm{cmd}
    &=
      \text{clamp}
      \left(
        \dfrac{P_0^\mathrm{ord}}{V_{\mathrm{safe},0}^\mathrm{meas}};\,
        0,\,
        I_{p,0}^{\max}
      \right)
\end{aligned}
```

A standard start requires $V_{T,0}>0$,
$V_{\mathrm{dip}} < V_{T,0} < V_{\mathrm{up}}$, and
$s_P^\mathrm{off} + s_P\omega_{g,0} \ne 0$. It also requires $Q_0^\mathrm{ref}$ within
$[Q^{\min},Q^{\max}]$, $P_0^\mathrm{ord}$ within
$[P^{\min},P^{\max}]$, both
$k_{\mathrm{base}}Q_0^\mathrm{gen}/V_{\mathrm{safe},0}^\mathrm{meas}$ and
$I_{q,0}^\mathrm{raw}$ within $[-I_{q,0}^{\max},I_{q,0}^{\max}]$, and
$P_0^\mathrm{ord}/V_{\mathrm{safe},0}^\mathrm{meas}$ within
$[0,I_{p,0}^{\max}]$. The priority-dependent current-circle radicand must be
nonnegative.

### Output Initialization

```math
\begin{aligned}
  P_e
    &\leftarrow V_{\mathrm{safe},0}^\mathrm{meas} I_{p,0}^\mathrm{cmd} \\
  Q^\mathrm{gen}
    &\leftarrow V_{\mathrm{safe},0}^\mathrm{meas} I_{q,0}^\mathrm{cmd} \\
  Q^\mathrm{ext}
    &\leftarrow Q_0^\mathrm{gen} \\
  \phi^\mathrm{ref}
    &\leftarrow
      \begin{cases}
        \tan^{-1}\!\left(Q_0^\mathrm{gen}/P_{e,0}\right) & P_{e,0} \ne 0 \\
        0 & P_{e,0} = 0
      \end{cases} \\
  P^\mathrm{ref}
    &\leftarrow
      \dfrac{P_{e,0}}
            {s_P^\mathrm{off} + s_P\omega_{g,0}}
\end{aligned}
```

REECA writes the resolved feedback and reference values to attached `pe`,
`qgen`, `qext`, `pfaref`, and `pref` signal inputs. If no signal is attached,
those values are used as constant inputs.

## Monitorable Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`iqcmd`         | [p.u.] | Reactive-current command output     | $I_q^\mathrm{cmd}$ (system base)
`ipcmd`         | [p.u.] | Active-current command output       | $I_p^\mathrm{cmd}$ (system base)
`vmeas`         | [p.u.] | Filtered terminal voltage           | $V^\mathrm{meas}$
`pmeas`         | [p.u.] | Filtered electrical power           | $P^\mathrm{meas}$ (component base)
