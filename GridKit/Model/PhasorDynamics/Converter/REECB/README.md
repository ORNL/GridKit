# **Renewable Energy Electrical Control Model (REECB)**

REECB is a WECC renewable electrical-control model for inverter-coupled
resources.

## Notes

- REECB is a control model only. It measures the terminal bus and publishes
  current commands; it injects no current into the network.
- When used with REPCA active-power control, connect REPCA `pext` to REECB `pref`.

## Block Diagram

![REECB electrical-control block diagram](../../../../../docs/Figures/PhasorDynamics/REECB/diagram.png)

Figure 1: REECB electrical-control model. Figure courtesy of the
[PowerWorld REEC_B model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20REEC_B.htm).

## Model Parameters

Symbol                              | Units    | JSON     | Description                                             | Typical Value | Note
------------------------------------|----------|----------|---------------------------------------------------------|---------------|------
$S^\mathrm{base}$                   | [MVA]    | `mva`    | REECB component power base                              | 100.0         | Required positive value; block name: `MVABase`
$s_{\mathrm{pf}}$                   | [binary] | `PfFlag` | Power-factor control flag                               | 0             | Block name: `PfFlag`; 1 = power-factor control, 0 = Q control
$s_V$                               | [binary] | `VFlag`  | Voltage-control mode flag                               | 0             | Block name: `VFlag`; 1 = Q control, 0 = voltage control
$s_Q$                               | [binary] | `QFlag`  | Reactive-power control flag                             | 0             | Block name: `QFlag`; 1 = voltage/Q control, 0 = constant pf or Q control
$s_{PQ}$                            | [binary] | `Pqflag` | P/Q priority flag for converter current limit           | 0             | Block name: `Pqflag`; 0 = Q priority, 1 = P priority
$T_{\mathrm{rv}}$                   | [sec]    | `Trv`    | Voltage-measurement filter time constant                | 0.02          | State 1; raised to the minimum-time floor
$T_{\mathrm{p}}$                    | [sec]    | `Tp`     | Electrical-power measurement filter time constant       | 0.0           | State 2; raised to the minimum-time floor
$V_0^\mathrm{ref}$                  | [p.u.]   | `Vref0`  | Outer-loop voltage reference                            | $V_{T,0}$     | Initialized from terminal voltage if omitted
$V_{\mathrm{dip}}$                  | [p.u.]   | `Vdip`   | Low-voltage threshold for the voltage-band gate         | 0.85          |
$V_{\mathrm{up}}$                   | [p.u.]   | `Vup`    | High-voltage threshold for the voltage-band gate        | 1.15          |
$D_1^\mathrm{db}$                   | [p.u.]   | `dbd1`   | Lower deadband threshold for voltage-error response     | 0.0           |
$D_2^\mathrm{db}$                   | [p.u.]   | `dbd2`   | Upper deadband threshold for voltage-error response     | 0.0           |
$K_{\mathrm{qv}}$                   | [p.u.]   | `kqv`    | Reactive-current injection gain outside the voltage band | 5.0          |
$I_{q,\mathrm{inj}}^{\min}$         | [p.u.]   | `Iql1`   | Minimum reactive-current injection limit                | -1.1          |
$I_{q,\mathrm{inj}}^{\max}$         | [p.u.]   | `Iqh1`   | Maximum reactive-current injection limit                | 1.1           |
$Q^{\max}$                          | [p.u.]   | `Qmax`   | Maximum reactive-power control limit                    | 0.436         |
$Q^{\min}$                          | [p.u.]   | `Qmin`   | Minimum reactive-power control limit                    | -0.436        |
$K_{\mathrm{qp}}$                   | [p.u.]   | `Kqp`    | Reactive-power control proportional gain                | 0.0           |
$K_{\mathrm{qi}}$                   | [p.u./s] | `Kqi`    | Reactive-power control integral gain                    | 0.1           |
$V^{\max}$                          | [p.u.]   | `Vmax`   | Maximum voltage-control limit                           | 1.1           |
$V^{\min}$                          | [p.u.]   | `Vmin`   | Minimum voltage-control limit                           | 0.9           |
$K_{\mathrm{vp}}$                   | [p.u.]   | `Kvp`    | Voltage-control proportional gain                       | 18.0          |
$K_{\mathrm{vi}}$                   | [p.u./s] | `Kvi`    | Voltage-control integral gain                           | 5.0           |
$T_{\mathrm{iq}}$                   | [sec]    | `Tiq`    | Reactive-current command lag time constant              | 0.02          | State 5; raised to the minimum-time floor
$T_{\mathrm{pord}}$                 | [sec]    | `Tpord`  | Active-power order filter time constant                 | 0.02          | State 6; raised to the minimum-time floor
$R_P^{\max}$                        | [p.u./s] | `dPmax`  | Positive active-power order ramp-rate limit             | 99.0          |
$R_P^{\min}$                        | [p.u./s] | `dPmin`  | Negative active-power order ramp-rate limit             | -99.0         |
$P^{\max}$                          | [p.u.]   | `Pmax`   | Maximum active-power order limit                        | 1.0           |
$P^{\min}$                          | [p.u.]   | `Pmin`   | Minimum active-power order limit                        | 0.0           |
$I^{\max}$                          | [p.u.]   | `Imax`   | Maximum total converter current                         | 1.3           |

### Parameter Validation

Invalid REECB parameter sets are rejected by the following checks:

```math
\begin{aligned}
  S^\mathrm{base} &> 0 \\
  s_{\mathrm{pf}}, s_V, s_Q, s_{PQ}
    &\in \{0,1\} \\
  T_{\mathrm{rv}}, T_{\mathrm{p}}, T_{\mathrm{iq}}, T_{\mathrm{pord}}
    &\ge 0 \\
  V_{\mathrm{dip}}
    &< V_{\mathrm{up}} \\
  D_1^\mathrm{db}
    &\le 0 \le D_2^\mathrm{db} \\
  I_{q,\mathrm{inj}}^{\min}
    &\le I_{q,\mathrm{inj}}^{\max} \\
  Q^{\min}
    &\le Q^{\max} \\
  V^{\min}
    &\le V^{\max} \\
  R_P^{\min}
    &< 0 < R_P^{\max} \\
  P^{\min}
    &\le P^{\max} \\
  I^{\max}
    &\ge 0
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!\left(T_x,\epsilon_T\right),
       \quad x\in\{\mathrm{rv},\mathrm{p},\mathrm{iq},\mathrm{pord}\} \\
  s_{\mathrm{pf}}^\mathrm{off}
    &= 1 - s_{\mathrm{pf}} \\
  s_V^\mathrm{off}
    &= 1 - s_V \\
  s_Q^\mathrm{off}
    &= 1 - s_Q \\
  k_\mathrm{base}
    &= \dfrac{S^\mathrm{sys}}{S^\mathrm{base}}
\end{aligned}
```

Multiplying by $k_\mathrm{base}$ converts system base to component base.

## Model Ports

Name     | Port   | Init    | Base        | Description
---------|--------|---------|-------------|------
`bus`    | Bus    | Known   | -           | Terminal-bus voltage
`pe`     | Input  | Unknown | System      | Electrical active-power feedback
`qgen`   | Input  | Unknown | System      | Reactive-power feedback
`qext`   | Input  | Unknown | System      | External reactive-power command
`pfaref` | Input  | Unknown | [rad]       | Power-factor angle reference
`pref`   | Input  | Unknown | System      | External active-power reference
`iqcmd`  | Output | Known   | System      | Reactive-current command
`ipcmd`  | Output | Known   | System      | Active-current command

`Known` ports are seeded before `initialize()` and preserved by it. `Unknown`
inputs are resolved during initialization and written to attached signal
storage, or retained as constant inputs when the port is unattached.

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
$s_{\mathrm{dip}}$                  | [-]      | Smooth voltage inside-band control gate | Approximately 1 inside the voltage band
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

Symbol                               | Units  | Init    | Description                       | Note
-------------------------------------|--------|---------|-----------------------------------|------
$V_{\mathrm{r}}$                     | [p.u.] | Known   | Terminal voltage, real component  | Bus input
$V_{\mathrm{i}}$                     | [p.u.] | Known   | Terminal voltage, imaginary component | Bus input
$P_e$                                | [p.u.] | Unknown | Electrical active-power feedback  | Signal port `pe`; system base
$Q^\mathrm{gen}$                     | [p.u.] | Unknown | Reactive-power feedback           | Signal port `qgen`; system base
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
      \left(k_\mathrm{base}P_e - P^\mathrm{meas}\right) \\
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

CommonMath defines the [`antiwindup`](../../../../CommonMath.md#antiwindup)
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
    + s_{\mathrm{pf}}^\mathrm{off}k_\mathrm{base}Q^\mathrm{ext} \\
  0 &=
    -e_Q
    + \text{clamp}
      \left(Q^\mathrm{ref};\, Q^{\min}, Q^{\max}\right)
    - k_\mathrm{base}Q^\mathrm{gen} \\
  0 &=
    -V_Q^\mathrm{PI}
    + \text{clamp}
      \left(K_{\mathrm{qp}}e_Q + x_Q^\mathrm{PI};\,
            V^{\min}, V^{\max}\right) \\
  0 &=
    -e_V^\mathrm{PI}
    + s_V V_Q^\mathrm{PI}
    + s_V^\mathrm{off}Q^\mathrm{ref}
    - V^\mathrm{meas} \\
  0 &=
    -f_P^\mathrm{ord}
    + \dfrac{1}{T_{\mathrm{pord}}}
      \left(k_\mathrm{base}P^\mathrm{ref} - P^\mathrm{ord}\right) \\
  0 &=
    -r_P^\mathrm{ord}
    + \text{clamp}
      \left(f_P^\mathrm{ord};\, R_P^{\min}, R_P^{\max}\right) \\
  0 &=
    -\left(I_q^\mathrm{circ}\right)^2
    + \left(I^{\max}\right)^2
    - s_{PQ}\left(k_\mathrm{base}I_p^\mathrm{cmd}\right)^2 \\
  0 &=
    -\left(I_p^\mathrm{circ}\right)^2
    + \left(I^{\max}\right)^2
    - \left(1-s_{PQ}\right)\left(k_\mathrm{base}I_q^\mathrm{cmd}\right)^2 \\
  0 &=
    -I_q^{\max}
    + \left(1-s_{PQ}\right)I^{\max}
    + s_{PQ}I_q^\mathrm{circ} \\
  0 &=
    -I_p^{\max}
    + s_{PQ}I^{\max}
    + \left(1-s_{PQ}\right)I_p^\mathrm{circ} \\
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
    -I_q^\mathrm{cmd}
    + \dfrac{1}{k_\mathrm{base}}\text{clamp}
      \left(I_q^\mathrm{raw};\,
            -I_q^{\max}, I_q^{\max}\right) \\
  0 &=
    -I_p^\mathrm{cmd}
    + \dfrac{1}{k_\mathrm{base}}\text{clamp}
      \left(
        \dfrac{P^\mathrm{ord}}{V_{\mathrm{safe}}^\mathrm{meas}};\,
        0,\,
        I_p^{\max}
      \right)
\end{aligned}
```

CommonMath defines helper targets and smooth approximations for
[max, clamp, deadband2, and inside](../../../../CommonMath.md#derived-functions).

## Initialization

### Input Initialization

The upstream source model seeds `ipcmd` and `iqcmd` before REECB initializes.
REECB snapshots them on component base first:

```math
\begin{aligned}
  V_{\mathrm{r}}, V_{\mathrm{i}}
    &\leftarrow \text{terminal-bus voltage} \\
  I_p^\mathrm{seed}
    &\leftarrow k_\mathrm{base}I_p^\mathrm{cmd} \\
  I_q^\mathrm{seed}
    &\leftarrow k_\mathrm{base}I_q^\mathrm{cmd}
\end{aligned}
```

Initialization never replaces the system-base values held in
$I_p^\mathrm{cmd}$ and $I_q^\mathrm{cmd}$.

### Internal Initialization

The residual limits with the smooth CommonMath
[`clamp`](../../../../CommonMath.md#clamp), so a steady state is seeded with
the limiter *input*, not its output. With initialization tolerance
$\epsilon_0=10^{-10}$, $\text{clamp}^{-1}(z;\ell,u)$ is the input producing
output $z$, and $u_0^\mathrm{aw}(a,f;\ell,u)$ the input holding an anti-windup
path stationary: $a$ when $|f|\le\epsilon_0$, else just past the limit $f$
drives toward. Both reject $z$ outside $[\ell,u]$.

Subscript $0$ denotes initial values; all internal derivatives start at zero:

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
    &= V_{\mathrm{safe},0}^\mathrm{meas}I_p^\mathrm{seed} \\
  k_\mathrm{base}Q_0^\mathrm{gen}
    &= V_{\mathrm{safe},0}^\mathrm{meas}I_q^\mathrm{seed} \\
  I_{q,0}^\mathrm{circ}
    &=
      \begin{cases}
        I^{\max} & s_{PQ}=0 \\
        \sqrt{(I^{\max})^2-(I_p^\mathrm{seed})^2}
          & s_{PQ}=1
      \end{cases} \\
  I_{p,0}^\mathrm{circ}
    &=
      \begin{cases}
        \sqrt{(I^{\max})^2-(I_q^\mathrm{seed})^2}
          & s_{PQ}=0 \\
        I^{\max} & s_{PQ}=1
      \end{cases} \\
  I_{q,0}^{\max}
    &= (1-s_{PQ})I^{\max}+s_{PQ}I_{q,0}^\mathrm{circ} \\
  I_{p,0}^{\max}
    &= s_{PQ}I^{\max}+(1-s_{PQ})I_{p,0}^\mathrm{circ}
\end{aligned}
```

```math
\begin{aligned}
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
  I_{q,0}^\mathrm{raw}
    &= \text{clamp}^{-1}
       \left(I_q^\mathrm{seed};\,
             -I_{q,0}^{\max},I_{q,0}^{\max}\right) \\
  I_{q,0}^\mathrm{control}
    &= I_{q,0}^\mathrm{raw}
       -(1-s_{\mathrm{dip},0})I_{q,0}^\mathrm{inj} \\
  u_{p,0}
    &= \text{clamp}^{-1}
       \left(I_p^\mathrm{seed};\,0,I_{p,0}^{\max}\right) \\
  P_0^\mathrm{ord}
    &= V_{\mathrm{safe},0}^\mathrm{meas}u_{p,0} \\
  f_{P,0}^\mathrm{ord}
    &= \text{clamp}^{-1}
       \left(0;\,R_P^{\min},R_P^{\max}\right) \\
  r_{P,0}^\mathrm{ord}
    &= \text{clamp}
       \left(f_{P,0}^\mathrm{ord};\,R_P^{\min},R_P^{\max}\right)
\end{aligned}
```

```math
\begin{aligned}
  Q_0^\mathrm{ref}
    &=
      \begin{cases}
        \text{clamp}^{-1}
        \left(k_\mathrm{base}Q_0^\mathrm{gen};\,
              Q^{\min},Q^{\max}\right)
          & s_Q=1\ \land\ s_V=1 \\
        V_0^\mathrm{meas}
          & s_Q=1\ \land\ s_V=0 \\
        V_{\mathrm{safe},0}^\mathrm{meas}I_{q,0}^\mathrm{control}
          & s_Q=0
      \end{cases} \\
  e_{Q,0}
    &=
      \text{clamp}
      \left(Q_0^\mathrm{ref};\, Q^{\min}, Q^{\max}\right)
      - k_\mathrm{base}Q_0^\mathrm{gen} \\
  Q_{V,0}
    &= \dfrac{Q_0^\mathrm{ref}}{V_{\mathrm{safe},0}^\mathrm{meas}} \\
  u_{Q,0}^\mathrm{PI}
    &=
      \begin{cases}
        \text{clamp}^{-1}
        \left(V_0^\mathrm{meas};\,V^{\min},V^{\max}\right)
          & s_Q=1\ \land\ s_V=1 \\
        u_0^\mathrm{aw}
        \left(Q_0^\mathrm{ref},K_{\mathrm{qi}}e_{Q,0};\,
              V^{\min},V^{\max}\right)
          & s_Q=1\ \land\ s_V=0 \\
        u_0^\mathrm{aw}
        \left(s_VV_0^\mathrm{meas}+s_V^\mathrm{off}Q_0^\mathrm{ref},\,
              K_{\mathrm{qi}}e_{Q,0};\,V^{\min},V^{\max}\right)
          & s_Q=0
      \end{cases} \\
  V_{Q,0}^\mathrm{PI}
    &=
      \text{clamp}
      \left(u_{Q,0}^\mathrm{PI};\, V^{\min}, V^{\max}\right) \\
  e_{V,0}^\mathrm{PI}
    &=
      s_V V_{Q,0}^\mathrm{PI}
      + s_V^\mathrm{off}Q_0^\mathrm{ref}
      - V_0^\mathrm{meas} \\
  x_{Q,0}^\mathrm{PI}
    &= u_{Q,0}^\mathrm{PI} - K_{\mathrm{qp}}e_{Q,0} \\
  u_{V,0}^\mathrm{PI}
    &=
      \begin{cases}
        \text{clamp}^{-1}
        \left(I_{q,0}^\mathrm{control};\,
              -I_{q,0}^{\max},I_{q,0}^{\max}\right)
          & s_Q=1 \\
        u_0^\mathrm{aw}
        \left(0,K_{\mathrm{vi}}e_{V,0}^\mathrm{PI};\,
              -I_{q,0}^{\max},I_{q,0}^{\max}\right)
          & s_Q=0
      \end{cases} \\
  I_{q,0}^\mathrm{base}
    &=
      \text{clamp}
      \left(u_{V,0}^\mathrm{PI};\,
            -I_{q,0}^{\max},
            I_{q,0}^{\max}\right) \\
  x_{V,0}^\mathrm{PI}
    &= u_{V,0}^\mathrm{PI} - K_{\mathrm{vp}}e_{V,0}^\mathrm{PI} \\
  0
    &= -I_{q,0}^\mathrm{raw}
       +s_Q I_{q,0}^\mathrm{base}
       +s_Q^\mathrm{off}Q_{V,0}
       +(1-s_{\mathrm{dip},0})I_{q,0}^\mathrm{inj}
\end{aligned}
```

The $s_Q=1$ path initializes the voltage PI output to reproduce
$I_{q,0}^\mathrm{control}$. The $s_Q=0$ path instead carries that target in
$Q_{V,0}$ and parks the otherwise inactive voltage PI path consistently.

Initialization rejects an operating point when any of the following holds:

- the bus voltage or either command seed is not finite;
- $I_p^\mathrm{seed}<0$;
- the command seeds leave the $I^{\max}$ circle, or a selected priority-circle
  radicand is less than $-\epsilon_0$;
- the physical active-power target
  $V_{\mathrm{safe},0}^\mathrm{meas}I_p^\mathrm{seed}$ lies outside
  $[P^{\min},P^{\max}]$ by more than $\epsilon_0$;
- a required current, ramp-rate, reactive-power, voltage, or controller output
  has no limiter input on its selected limits; or
- $s_{\mathrm{pf}}=1$, $|P_0^\mathrm{meas}|\le\epsilon_0$, and
  $|Q_0^\mathrm{ref}|>\epsilon_0$.

Every check resolves before any storage is written, so a rejected
initialization leaves state, command nodes, and external signals unchanged.

### Output Initialization

```math
\begin{aligned}
  P_{e,0}
    &\leftarrow V_{\mathrm{safe},0}^\mathrm{meas}I_p^\mathrm{cmd} \\
  Q_0^\mathrm{gen}
    &\leftarrow V_{\mathrm{safe},0}^\mathrm{meas}I_q^\mathrm{cmd} \\
  Q_0^\mathrm{ext}
    &\leftarrow \dfrac{Q_0^\mathrm{ref}}{k_\mathrm{base}} \\
  \phi_0^\mathrm{ref}
    &\leftarrow
      \begin{cases}
        \tan^{-1}\!\left(Q_0^\mathrm{ref}/P_0^\mathrm{meas}\right)
          & s_{\mathrm{pf}}=1
            \ \land\ |P_0^\mathrm{meas}|>\epsilon_0 \\
        0
          & s_{\mathrm{pf}}=0
            \ \lor\
            \left(|P_0^\mathrm{meas}|\le\epsilon_0
            \ \land\ |Q_0^\mathrm{ref}|\le\epsilon_0\right)
      \end{cases} \\
  P_0^\mathrm{ref}
    &\leftarrow
      \dfrac{P_0^\mathrm{ord}
             +T_{\mathrm{pord}}f_{P,0}^\mathrm{ord}}
            {k_\mathrm{base}}
\end{aligned}
```

These expressions are on system base. $Q_0^\mathrm{ext}$ is published even when
$s_{\mathrm{pf}}=1$, though the power-factor path does not consume it.

## Monitorable Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`iqcmd`         | [p.u.] | Reactive-current command output     | $I_q^\mathrm{cmd}$ (system base)
`ipcmd`         | [p.u.] | Active-current command output       | $I_p^\mathrm{cmd}$ (system base)
`vmeas`         | [p.u.] | Filtered terminal voltage           | $V^\mathrm{meas}$
`pmeas`         | [p.u.] | Filtered electrical power           | $P^\mathrm{meas}$ (component base)
