# GGOV1

GGOV1 is a general governor-turbine model with electrical-power measurement,
speed/load reference selection, proportional/integral/derivative governor
control, load limiting, acceleration limiting, temperature limiting, actuator
rate limits, turbine lag/lead dynamics, and optional diesel damping.

> [!WARNING]
> The damping, fuel-flow, PID, limiter, temperature, and transport-delay paths
> remain incomplete or unverified. The initialization below does not yet
> establish a consistent steady state for the smooth selector paths.

## Notes

- Internal control, valve-stroke, and turbine-power quantities are on the
  GGOV1 component base unless otherwise stated.
- The dashed speed deadband block and `Db` source field are only for GGOV1D.
  GGOV1 uses the speed input directly.
- The source diagram notes that `Rup` and `Rdown` inputs are not implemented in
  Simulator; the equations below do not use those source fields.

## Block Diagram

![](../../../../../docs/Figures/PhasorDynamics/GGOV1_diagram.png)

Figure 1: Governor GGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol              | Units     | JSON      | Description                                     | Typical Value | Note
--------------------|-----------|-----------|-------------------------------------------------|---------------|----------------------------------------------------------------------------------------------------
$P^{\mathrm{rate}}$ | [MW]      | `Trate`   | Optional turbine-rating power base              | 0.0           | `Trate > 0` defines the governor base
$I_R$               | [integer] | `Rselect` | Droop feedback selector                         | 1             | selects speed, electrical power, or valve feedback
$s_\mathrm{flag}$   | [binary]  | `Flag`    | Turbine-speed multiplier selector               | 1             | 1 uses $1+\omega$, 0 uses 1.0
$R$                 | [p.u.]    | `R`       | Permanent droop                                 | 0.05          | Source label: `r`
$T_\mathrm{pelec}$  | [s]       | `Tpelec`  | Electrical-power measurement time constant      | 0.0           | State 1 in Fig. 1
$e^{\max}$          | [p.u.]    | `Maxerr`  | Upper limit on the signed governor error        | 1.0           | Source label: `maxerr`
$e^{\min}$          | [p.u.]    | `Minerr`  | Lower limit on the signed governor error        | -1.0          | Source label: `minerr`
$K_\mathrm{pgov}$   | [p.u.]    | `Kpgov`   | Governor proportional gain                      | 10.0          |
$K_\mathrm{igov}$   | [p.u./s]  | `Kigov`   | Governor integral gain                          | 1.0           | State 3
$K_\mathrm{dgov}$   | [p.u.]    | `Kdgov`   | Governor differential gain                      | 0.0           | State 2
$T_\mathrm{dgov}$   | [s]       | `Tdgov`   | Governor differential time constant             | 0.0           |
$V^{\max}$          | [p.u.]    | `Vmax`    | Maximum governor output before actuator         | 1.0           | Governor response limits may adjust this value
$V^{\min}$          | [p.u.]    | `Vmin`    | Minimum governor output before actuator         | 0.0           | Governor response limits may adjust this value
$T_\mathrm{act}$    | [s]       | `Tact`    | Turbine actuator time constant                  | 0.1           | State 4 in Fig. 1
$R_\mathrm{open}$   | [p.u./s]  | `Ropen`   | Maximum actuator opening rate                   | 1.0           |
$R_\mathrm{close}$  | [p.u./s]  | `Rclose`  | Maximum actuator closing rate                   | -1.0          |
$K_\mathrm{turb}$   | [p.u.]    | `Kturb`   | Turbine gain                                    | 1.0           |
$W_\mathrm{fnl}$    | [p.u.]    | `Wfnl`    | No-load fuel flow                               | 0.0           |
$T_B$               | [s]       | `Tb`      | Turbine lead–lag denominator time constant      | 0.0           | State 5 in Fig. 1
$T_C$               | [s]       | `Tc`      | Turbine lead–lag numerator time constant        | 0.0           |
$T_\mathrm{eng}$    | [s]       | `Teng`    | Engine transport lag                            | 0.0           | Source label: `e^{-sTeng}`; source transport delay is not represented as a differential state below
$T_\mathrm{fload}$  | [s]       | `Tfload`  | Load-limiter lag time constant                  | 0.0           | State 6 in Fig. 1
$K_\mathrm{pload}$  | [p.u.]    | `Kpload`  | Load-limiter proportional gain                  | 0.0           | note path changes when zero
$K_\mathrm{iload}$  | [p.u./s]  | `Kiload`  | Load-limiter integral gain                      | 0.0           | State 7 in Fig. 1
$L_\mathrm{dref}$   | [p.u.]    | `Ldref`   | Load reference                                  | 1.0           |
$D_\mathrm{m}$               | [p.u.]    | `Dm`      | Diesel damping gain                             | 0.0           | sign-dependent speed term in Fig. 1
$K_\mathrm{imw}$    | [p.u./s]  | `Kimw`    | Supervisory load-control integral gain          | 0.0           | State 8 in Fig. 1
$A_\mathrm{set}$    | [p.u.]    | `Aset`    | Acceleration-control reference                  | 0.0           | Source label: `aset`
$K_A$               | [p.u.]    | `Ka`      | Acceleration-control gain                       | 0.0           | Source label: `KA`
$T_A$               | [s]       | `Ta`      | Acceleration-control time constant              | 0.0           | State 9 in Fig. 1
$T_\mathrm{sa}$     | [s]       | `Tsa`     | Temperature-detection numerator time constant   | 0.0           | State 10 in Fig. 1
$T_\mathrm{sb}$     | [s]       | `Tsb`     | Temperature-detection denominator time constant | 0.0           | State 10 in Fig. 1
$R_\mathrm{up}$     | [p.u./s]  | `Rup`     | Source upward ramp input                        | 0.0           | Source note says not implemented in Simulator
$R_\mathrm{down}$   | [p.u./s]  | `Rdown`   | Source downward ramp input                      | 0.0           | Source note says not implemented in Simulator

### Parameter Validation

A valid GGOV1 parameter set must satisfy the following conditions:

```math
\begin{aligned}
  &P^{\mathrm{rate}}\ge 0,\quad R\gt 0,\quad I_R\in\{-2,-1,1\},\quad s_\mathrm{flag}\in\{0,1\} \\
  &T_\mathrm{pelec},T_\mathrm{dgov},T_B,T_C,T_\mathrm{eng},T_\mathrm{fload},T_A,T_\mathrm{sa},T_\mathrm{sb}\ge 0,\quad T_\mathrm{act}\gt 0 \\
  &T_B \gt 0\quad\text{or}\quad(T_B = 0\ \text{and}\ T_C = 0) \\
  &e^{\min}\le e^{\max},\quad V^{\min}\le V^{\max},\quad R_\mathrm{close}\lt 0\lt R_\mathrm{open} \\
  &K_\mathrm{turb}\gt 0
\end{aligned}
```

### Model Derived Parameters

The component base and flag complements are:

```math
\begin{aligned}
  S^{\mathrm{base}}
    &=
    \begin{cases}
      P^{\mathrm{rate}} & P^{\mathrm{rate}} \gt 0 \\
      S^{\mathrm{sys}} & \text{otherwise}
    \end{cases} \\
  s_\mathrm{flag}^{\mathrm{off}} &= 1 - s_\mathrm{flag}
\end{aligned}
```

$S^{\mathrm{sys}}$ is the system power base.

## Model Ports

Name       | Port   | Init | Description
-----------|--------|------|------------
`pref`     | Input  | TBD  | Governor reference
`paux`     | Input  | TBD  | Auxiliary power input
`pmwset`   | Input  | TBD  | Supervisory MW setpoint
`pelec`    | Input  | TBD  | Electrical active power
`ldref`    | Input  | TBD  | Load reference
`speed`    | Input  | TBD  | Machine speed deviation
`pmech`    | Output | TBD  | Mechanical-power output

## Model Variables

### Internal Variables

#### Differential

Symbol                            | Units  | Description                          | Note
----------------------------------|--------|--------------------------------------|--------------------------------------------------
$P_\mathrm{elec}^{\mathrm{meas}}$ | [p.u.] | Measured electrical power            | State 1 in Fig. 1; Source label: `Pelec Measured`
$x_D$                             | [p.u.] | Governor differential control state  | State 2 in Fig. 1
$x_I$                             | [p.u.] | Governor integral control state      | State 3 in Fig. 1
$x_\mathrm{act}$                  | [p.u.] | Turbine actuator or valve stroke     | State 4 in Fig. 1
$x_\mathrm{turb}$                 | [p.u.] | Turbine lead–lag state               | State 5 in Fig. 1; Source label: `Turbine LL`
$x_\mathrm{load}$                 | [p.u.] | Turbine load-limiter lag state       | State 6 in Fig. 1
$x_\mathrm{ldint}$                | [p.u.] | Turbine load integral-control state  | State 7 in Fig. 1
$x_\mathrm{mw}$                   | [p.u.] | Supervisory load-control state       | State 8 in Fig. 1
$x_\mathrm{acc}$                  | [p.u.] | Acceleration-control state           | State 9 in Fig. 1
$x_\mathrm{temp}$                 | [p.u.] | Temperature-detection lead–lag state | State 10 in Fig. 1

#### Algebraic

Symbol             | Units    | Description                               | Note
-------------------|----------|-------------------------------------------|----------------------------------------------------
$P_\mathrm{mwref}$ | [p.u.]   | Supervisory load-control reference        | From $P_\mathrm{mwset}-P_\mathrm{elec}$
$y_R$              | [p.u.]   | Selected droop feedback                   | Controlled by `Rselect`
$e_G$              | [p.u.]   | Limited governor error                    | After $e^{\min}$ and $e^{\max}$
$f_\mathrm{pid}$   | [p.u.]   | Governor PID output                       | Forms `fsrn`
$f_\mathrm{srn}$   | [p.u.]   | Normal governor fuel/stroke request       | Low-value select input
$f_\mathrm{sra}$   | [p.u.]   | Acceleration-control request              | Low-value select input
$f_\mathrm{srt}$   | [p.u.]   | Temperature/load request                  | Low-value select input
$f_\mathrm{srl}$   | [p.u.]   | Acceleration/temperature low-value select | Lesser of $f_\mathrm{sra}$ and $f_\mathrm{srt}$
$f_\mathrm{sr}$    | [p.u.]   | Low-value select output                   | Limited by $V^{\min}$ and $V^{\max}$
$r_\mathrm{act}$   | [p.u./s] | Actuator rate-limited derivative          | Limited by $R_\mathrm{close}$ and $R_\mathrm{open}$
$P_\mathrm{turb}$  | [p.u.]   | Turbine power before damping              | After turbine lead–lag and transport lag
$P_\mathrm{damp}$  | [p.u.]   | Damping power term                        | Source label: `Dm`
$P_\mathrm{m}$              | [p.u.]   | Mechanical-power output                   | Source label: `Pmech`

### External Variables

#### Differential

None.

#### Algebraic

Symbol                           | Units    | Description                           | Note
-------------------------------- | -------- | ------------------------------------- | ------------------------------------------------
$P_\mathrm{ref}$                 | [p.u.]   | Governor reference                    | Source label: `Pref`
$P_\mathrm{aux}$                 | [p.u.]   | Auxiliary power input                 | Source label: `Paux`; optional, defaults to zero
$P_\mathrm{mwset}$               | [p.u.]   | Supervisory MW setpoint               | Source label: `Pmwset`
$P_\mathrm{elec}$                | [p.u.]   | Electrical active power               | Source label: `Pelec`
$L_\mathrm{dref}$                | [p.u.]   | Load reference input                  | Source label: `Ldref`
$\omega$                         | [p.u.]   | Machine speed deviation               | Source label: `Speed`

## Model Equations

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup), [`clamp`](../../../../CommonMath.md#clamp), [`min`](../../../../CommonMath.md#minimum).

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -T_\mathrm{pelec}\dot P_\mathrm{elec}^{\mathrm{meas}} - P_\mathrm{elec}^{\mathrm{meas}} + P_\mathrm{elec} \\
  0 &= -T_\mathrm{dgov}\dot x_D - x_D + e_G \\
  0 &=
    -\dot x_I
    + \text{antiwindup}(
        f_\mathrm{pid},
        K_\mathrm{igov}e_G;
        V^{\min},
        V^{\max}
      ) \\
  0 &= -\dot x_\mathrm{act} + r_\mathrm{act} \\
  0 &= -T_B\dot x_\mathrm{turb} - x_\mathrm{turb} + x_\mathrm{act} \\
  0 &= -T_\mathrm{fload}\dot x_\mathrm{load} - x_\mathrm{load} + f_\mathrm{srt} \\
  0 &= -\dot x_\mathrm{ldint} + K_\mathrm{iload}(L_\mathrm{dref}-x_\mathrm{load}) \\
  0 &= -\dot x_\mathrm{mw} + K_\mathrm{imw}(P_\mathrm{mwset}-P_\mathrm{elec}) \\
  0 &= -T_A\dot x_\mathrm{acc} - x_\mathrm{acc} + \omega \\
  0 &= -T_\mathrm{sb}\dot x_\mathrm{temp} - x_\mathrm{temp} + f_\mathrm{sr}
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -P_\mathrm{mwref} + x_\mathrm{mw} + P_\mathrm{ref} + P_\mathrm{aux} \\
  0 &= -R y_R
       + \begin{cases}
          \omega & I_R = 1 \\
          P_\mathrm{elec}^{\mathrm{meas}} & I_R = -1 \\
          x_\mathrm{act} & I_R = -2
       \end{cases} \\
  0 &= -e_G + \text{clamp}(P_\mathrm{mwref} - y_R;\ e^{\min},\ e^{\max}) \\
  0 &= -f_\mathrm{pid} + K_\mathrm{pgov}e_G + K_\mathrm{dgov}(e_G - x_D) + x_I \\
  0 &= -f_\mathrm{srn} + \text{clamp}(f_\mathrm{pid}; V^{\min}, V^{\max}) \\
  0 &= -f_\mathrm{sra} + \text{clamp}(A_\mathrm{set} - K_A x_\mathrm{acc}; V^{\min}, V^{\max}) \\
  0 &= -f_\mathrm{srt} + \text{clamp}\left(\dfrac{L_\mathrm{dref} + P_\mathrm{aux}}{K_\mathrm{turb}} + W_\mathrm{fnl} + x_\mathrm{ldint}; V^{\min}, V^{\max}\right) \\
  0 &= -f_\mathrm{srl} + \text{min}(f_\mathrm{sra}, f_\mathrm{srt}) \\
  0 &= -f_\mathrm{sr} + \text{min}(f_\mathrm{srn}, f_\mathrm{srl}) \\
  0 &= -r_\mathrm{act} + \text{clamp}\left(\dfrac{f_\mathrm{sr}-x_\mathrm{act}}{T_\mathrm{act}}; R_\mathrm{close}, R_\mathrm{open}\right) \\
  0 &= -P_\mathrm{turb}
       + K_\mathrm{turb}
       \begin{cases}
         x_\mathrm{act} - W_\mathrm{fnl} & T_B = T_C = 0 \\
         x_\mathrm{turb} + \dfrac{T_C}{T_B}(x_\mathrm{act}-x_\mathrm{turb}) - W_\mathrm{fnl} & T_B \gt 0
       \end{cases} \\
  0 &= -P_\mathrm{damp}
       + D_\mathrm{m}
       \begin{cases}
          \omega & D_\mathrm{m} \ge 0 \\
          (1+\omega)^{D_\mathrm{m}} & D_\mathrm{m} \lt 0
       \end{cases} \\
  0 &= -P_\mathrm{m} + P_\mathrm{turb} + P_\mathrm{damp}
\end{aligned}
```

When $T_B=T_C=0$, the turbine lead–lag block is bypassed before the turbine
gain and no-load fuel-flow calculation.

### External Equations

None.

## Initialization

For a power-flow start, require $\omega=0$ and initialize all internal derivatives to zero:

```math
\begin{aligned}
  P_{\mathrm{aux}} &\leftarrow 0,\quad \text{if omitted} \\
  P_{\mathrm{elec}}^{\mathrm{meas}} &\leftarrow P_{\mathrm{elec}} \\
  x_{\mathrm{acc}} &\leftarrow 0 \\
  P_{\mathrm{damp}} &\leftarrow 0
\end{aligned}
```

Given initialized machine mechanical power, solve the actuator and turbine path:

```math
\begin{aligned}
  P_{\mathrm{turb}} &\leftarrow P_{\mathrm{m}} - P_{\mathrm{damp}} \\
  x_{\mathrm{act}} &\leftarrow W_\mathrm{fnl} + \dfrac{P_{\mathrm{turb}}}{K_\mathrm{turb}} \\
  x_{\mathrm{turb}} &\leftarrow x_{\mathrm{act}} \\
  f_{\mathrm{sr}} &\leftarrow x_{\mathrm{act}} \\
  f_{\mathrm{srn}},f_{\mathrm{sra}},f_{\mathrm{srt}},f_{\mathrm{srl}} &\leftarrow f_{\mathrm{sr}}
\end{aligned}
```

Initialize the limiter and control states:

```math
\begin{aligned}
  x_{\mathrm{load}} &\leftarrow f_{\mathrm{srt}} \\
  x_{\mathrm{ldint}}
    &\leftarrow f_{\mathrm{srt}}
       - \dfrac{L_{\mathrm{dref}}+P_{\mathrm{aux}}}{K_\mathrm{turb}}
       - W_\mathrm{fnl} \\
  x_{\mathrm{mw}} &\leftarrow 0 \\
  y_R &\leftarrow \dfrac{1}{R}
       \begin{cases}
          \omega & I_R = 1 \\
          P_{\mathrm{elec}}^{\mathrm{meas}} & I_R = -1 \\
          x_{\mathrm{act}} & I_R = -2
       \end{cases} \\
  P_{\mathrm{mwref}} &\leftarrow f_{\mathrm{pid}} + y_R \\
  P_{\mathrm{ref}} &\leftarrow P_{\mathrm{mwref}} - P_{\mathrm{aux}} - x_{\mathrm{mw}}
\end{aligned}
```

These assignments assume inactive low-value select alternatives,
inactive actuator rate limits, $V^{\min}\le f_{\mathrm{sr}}\le V^{\max}$,
and $K_\mathrm{turb}\ne 0$.

## Monitors

Monitor         | Units  | Description                          | Note
----------------|--------|--------------------------------------|------------------------
`pmech`         | [p.u.] | Mechanical-power output              | $P_\mathrm{m}$
`pelec_meas`    | [p.u.] | Measured electrical power            | State 1
`xd`            | [p.u.] | Governor differential-control state  | State 2
`xi`            | [p.u.] | Governor integral-control state      | State 3
`valve`         | [p.u.] | Turbine actuator or valve stroke     | State 4
`turbine_ll`    | [p.u.] | Turbine lead–lag state               | State 5
`load_limiter`  | [p.u.] | Turbine load-limiter state           | State 6
`load_int`      | [p.u.] | Turbine load integral-control state  | State 7
`mw_control`    | [p.u.] | Supervisory load-control state       | State 8
`accel_control` | [p.u.] | Acceleration-control state           | State 9
`temp_ll`       | [p.u.] | Temperature-detection lead–lag state | State 10
`fsrn`          | [p.u.] | Normal governor request              | Low-value select input
`fsra`          | [p.u.] | Acceleration-control request         | Low-value select input
`fsrt`          | [p.u.] | Temperature/load request             | Low-value select input
`fsr`           | [p.u.] | Selected governor request            | Low-value select output
