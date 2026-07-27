# **Renewable Energy Plant Control Model (REPCA)**

REPCA is a WECC renewable energy plant control model for inverter-coupled
resources.

## Notes

- REPCA is a control model only. It measures the regulated bus and branch
  signals, publishes plant power commands, and injects no current into the
  network.
- Branch-current, branch-power, power-reference, and command signals are on
  system base. Internal power states and controller limits are on component
  base.

## Block Diagram

![REPCA plant-control block diagram](../../../../../docs/Figures/PhasorDynamics/REPCA/diagram.png)

Figure 1: REPCA plant-control model. Figure courtesy of
[PowerWorld](https://www.powerworld.com/WebHelp/).

## Model Parameters

Symbol                              | Units    | JSON        | Description                                             | Typical Value | Note
------------------------------------|----------|-------------|---------------------------------------------------------|---------------|------
$S^\text{base}$                     | [MVA]    | `mva`       | REPCA component power base                              | 100.0         |
$s_\text{comp}$                     | [binary] | `VcompFlag` | Voltage-compensation mode flag                          | 1             | 1 = line-drop compensation, 0 = reactive droop
$s_\text{ref}$                      | [binary] | `RefFlag`   | Reactive-loop reference flag                            | 1             | 1 = voltage control, 0 = reactive-power control
$s_\text{freq}$                     | [binary] | `Freqflag`  | Active-power control output flag                        | 0             | 1 = command enabled, 0 = zero output
$T_\text{fltr}$                     | [sec]    | `Tfltr`     | Voltage and reactive-power measurement filter time constant | 0.05      |
$V^\text{frz}$                      | [p.u.]   | `Vfrz`      | Regulated-voltage threshold below which the reactive-power PI state freezes | 0.7 |
$R_c$                               | [p.u.]   | `Rc`        | Line-drop compensation resistance                       | 0.0           |
$X_c$                               | [p.u.]   | `Xc`        | Line-drop compensation reactance                        | 0.0           |
$K_c$                               | [p.u.]   | `Kc`        | Reactive-current compensation gain                      | 1.0           |
$D_\text{bd1}$                      | [p.u.]   | `dbdlow`    | Lower reactive-loop deadband threshold                  | 0.0           |
$D_\text{bd2}$                      | [p.u.]   | `dbdupper`  | Upper reactive-loop deadband threshold                  | 0.0           |
$e^{\max}$                          | [p.u.]   | `emax`      | Maximum reactive-loop error limit                       | 1.0           |
$e^{\min}$                          | [p.u.]   | `emin`      | Minimum reactive-loop error limit                       | -1.0          |
$K_\text{p}$                        | [p.u.]   | `Kp`        | Reactive-power controller proportional gain             | 10.0          |
$K_\text{i}$                        | [p.u./s] | `Ki`        | Reactive-power controller integral gain                 | 10.0          |
$Q^{\max}$                          | [p.u.]   | `Qmax`      | Maximum reactive-power command                          | 1.0           |
$Q^{\min}$                          | [p.u.]   | `Qmin`      | Minimum reactive-power command                          | -1.0          |
$T_\text{ft}$                       | [sec]    | `Tft`       | Reactive-command lead time constant                     | 0.0           |
$T_\text{fv}$                       | [sec]    | `Tfv`       | Reactive-command lag time constant                      | 3.0           |
$T_\text{p}$                        | [sec]    | `Tp`        | Active-power measurement filter time constant           | 0.0           |
$D_\text{bd1}^f$                    | [p.u.]   | `fdbd1`     | Lower frequency-error deadband threshold                | 0.0           |
$D_\text{bd2}^f$                    | [p.u.]   | `fdbd2`     | Upper frequency-error deadband threshold                | 0.0           |
$D_\text{dn}$                       | [p.u.]   | `Ddn`       | Down-frequency droop response gain                      | 20.0          |
$D_\text{up}$                       | [p.u.]   | `Dup`       | Up-frequency droop response gain                        | 0.0           |
$e_P^{\max}$                        | [p.u.]   | `femax`     | Maximum active-power error limit                        | 1.0           |
$e_P^{\min}$                        | [p.u.]   | `femin`     | Minimum active-power error limit                        | -1.0          |
$K_\text{pg}$                       | [p.u.]   | `Kpg`       | Active-power controller proportional gain               | 10.0          |
$K_\text{ig}$                       | [p.u./s] | `Kig`       | Active-power controller integral gain                   | 10.0          |
$P^{\max}$                          | [p.u.]   | `Pmax`      | Maximum active-power command                            | 2.0           |
$P^{\min}$                          | [p.u.]   | `Pmin`      | Minimum active-power command                            | 0.0           |
$T_\text{lag}$                      | [sec]    | `Tlag`      | Active-power command lag time constant                  | 3.0           |

### Parameter Validation

Invalid REPCA parameter sets are rejected by the following checks:

```math
\begin{aligned}
  S^\text{base}
    &> 0 \\
  s_\text{comp}, s_\text{ref}, s_\text{freq}
    &\in \{0,1\} \\
  T_\text{fv}
    &> 0 \\
  D_\text{bd1}
    &\le 0 \le D_\text{bd2} \\
  e^{\min}
    &\le 0 \le e^{\max} \\
  Q^{\min}
    &\le Q^{\max} \\
  D_\text{bd1}^f
    &\le 0 \le D_\text{bd2}^f \\
  e_P^{\min}
    &\le 0 \le e_P^{\max} \\
  P^{\min}
    &\le P^{\max}
\end{aligned}
```

`VcompFlag`, `RefFlag`, and `Freqflag` accept Boolean values or integer/real
values equal to 0 or 1.

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \text{s}$. Each time constant below $\epsilon_T$ is
raised to that floor and logged as a warning:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!\left(T_x,\epsilon_T\right),
       \quad x\in\{\text{fltr},\text{p},\text{lag}\} \\
  s_\text{comp}^\text{off} &= 1 - s_\text{comp} \\
  s_\text{ref}^\text{off} &= 1 - s_\text{ref} \\
  k_\text{base} &= \dfrac{S^\text{sys}}{S^\text{base}}
\end{aligned}
```

Multiplying by $k_\text{base}$ converts system base to component base.

## Model Ports

Name        | Port   | Init    | Base   | Description
------------|--------|---------|--------|------
`bus`       | Bus    | Known   | -      | Regulated-bus voltage
`ibranchr`  | Input  | Known   | System | Branch-current real component input
`ibranchi`  | Input  | Known   | System | Branch-current imaginary component input
`pbranch`   | Input  | Known   | System | Branch active-power input
`qbranch`   | Input  | Known   | System | Branch reactive-power input
`freq`      | Input  | Known   | -      | Frequency input
`freqref`   | Input  | Unknown | -      | Frequency reference
`vref`      | Input  | Unknown | -      | Voltage-control reference
`qref`      | Input  | Unknown | System | Reactive-power reference
`pplantref` | Input  | Unknown | System | Plant active-power reference
`qext`      | Output | Known   | System | Reactive-power command output
`pext`      | Output | Known   | System | Active-power command output

The known bus and input ports are seeded before `initialize()` and preserved by
it. `qext` is also seeded and preserved. `pext` is preserved when
$s_\text{freq}=1$ and set to zero when $s_\text{freq}=0$. Unknown reference
inputs are resolved during initialization and written to attached signal
storage, or retained as constant inputs when the port is unattached.

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                         | Note
------------------------|--------|-------------------------------------|------
$V^\text{meas}$         | [p.u.] | Filtered regulated voltage          | State 1 in Fig. 1
$Q^\text{meas}$         | [p.u.] | Filtered reactive-power signal      | State 2 in Fig. 1; component base
$x_Q^\text{PI}$         | [p.u.] | Reactive-power PI controller state  | State 3 in Fig. 1; component base
$x_Q^\text{lag}$        | [p.u.] | Reactive-command lead-lag state     | State 4 in Fig. 1; component base
$P^\text{meas}$         | [p.u.] | Filtered active-power signal        | State 5 in Fig. 1; component base
$x_P^\text{PI}$         | [p.u.] | Active-power PI controller state    | State 6 in Fig. 1; component base
$P^\text{ref}$          | [p.u.] | Active-power command lag state      | State 7 in Fig. 1; component base

#### Algebraic

Symbol                    | Units  | Description                         | Note
--------------------------|--------|-------------------------------------|------
$V$                       | [p.u.] | Regulated-bus voltage magnitude     |
$V^\text{ldc}$            | [p.u.] | Line-drop compensated voltage magnitude |
$V^\text{droop}$          | [p.u.] | Reactive-droop-compensated voltage |
$V^\text{ctrl}$           | [p.u.] | Selected voltage-measurement input  |
$s_\text{frz}$            | [-]    | Smooth reactive-power PI voltage-enable gate |
$e_\text{RQ}$             | [p.u.] | Selected reactive-loop error        |
$e_\text{RQ}^\text{db}$   | [p.u.] | Deadbanded reactive-loop error      |
$e_\text{RQ}^\text{lim}$  | [p.u.] | Limited reactive-loop error         |
$Q^\text{PI}$             | [p.u.] | Reactive-power PI output            | Component base
$Q^\text{ext}$            | [p.u.] | Reactive-power command output       | System base
$e_f$                     | [p.u.] | Frequency error after deadband      |
$e_P$                     | [p.u.] | Active-power control error          | Component base
$e_P^\text{lim}$          | [p.u.] | Limited active-power control error  | Component base
$P^\text{PI}$             | [p.u.] | Active-power PI output              | Component base
$P^\text{ext}$            | [p.u.] | Active-power command output         | System base

### External Variables

#### Differential
None.

#### Algebraic

Symbol                         | Units  | Init    | Description                       | Note
-------------------------------|--------|---------|-----------------------------------|------
$V_\text{r}$                   | [p.u.] | Known   | Regulated-bus voltage, real component | Bus input
$V_\text{i}$                   | [p.u.] | Known   | Regulated-bus voltage, imaginary component | Bus input
$I_\text{r}$                   | [p.u.] | Known   | Branch-current real component     | Signal port `ibranchr`; system base
$I_\text{i}$                   | [p.u.] | Known   | Branch-current imaginary component | Signal port `ibranchi`; system base
$P^\text{br}$                  | [p.u.] | Known   | Branch active power               | Signal port `pbranch`; system base
$Q^\text{br}$                  | [p.u.] | Known   | Branch reactive power             | Signal port `qbranch`; system base
$f$                            | [p.u.] | Known   | Frequency input                   | Signal port `freq`
$f^\text{ref}$                 | [p.u.] | Unknown | Frequency reference               | Optional signal port `freqref`
$V^\text{ref}$                 | [p.u.] | Unknown | Voltage-control reference         | Optional signal port `vref`
$Q^\text{ref}$                 | [p.u.] | Unknown | Reactive-power reference          | Optional signal port `qref`; system base
$P_\text{plant}^\text{ref}$    | [p.u.] | Unknown | Plant active-power reference      | Optional signal port `pplantref`; system base

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{V}^\text{meas}
    + \dfrac{1}{T_\text{fltr}}
      \left(V^\text{ctrl} - V^\text{meas}\right) \\
  0 &=
    -\dot{Q}^\text{meas}
    + \dfrac{1}{T_\text{fltr}}
      \left(k_\text{base}Q^\text{br} - Q^\text{meas}\right) \\
  0 &=
    -\dot{x}_Q^\text{PI}
    + s_\text{frz}\,
      \text{antiwindup}\!\left(
        Q^\text{PI},
        K_\text{i}e_\text{RQ}^\text{lim};\,
        Q^{\min},
        Q^{\max}
      \right) \\
  0 &=
    -\dot{x}_Q^\text{lag}
    + \dfrac{1}{T_\text{fv}}
      \left(Q^\text{PI} - x_Q^\text{lag}\right) \\
  0 &=
    -\dot{P}^\text{meas}
    + \dfrac{1}{T_\text{p}}
      \left(k_\text{base}P^\text{br} - P^\text{meas}\right) \\
  0 &=
    -\dot{x}_P^\text{PI}
    + \text{antiwindup}\!\left(
      P^\text{PI},
      K_\text{ig}e_P^\text{lim};\,
      P^{\min},
      P^{\max}
    \right) \\
  0 &=
    -\dot{P}^\text{ref}
    + \dfrac{1}{T_\text{lag}}
      \left(P^\text{PI} - P^\text{ref}\right)
\end{aligned}
```

CommonMath defines the [`antiwindup`](../../../../CommonMath.md#antiwindup)
target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -V^2
       + V_\text{r}^2
       + V_\text{i}^2 \\
  0 &=
    -\left(V^\text{ldc}\right)^2
    + \left(V_\text{r} - R_c I_\text{r} + X_c I_\text{i}\right)^2
    + \left(V_\text{i} - R_c I_\text{i} - X_c I_\text{r}\right)^2 \\
  0 &= -V^\text{droop} + V + K_c k_\text{base}Q^\text{br} \\
  0 &= -V^\text{ctrl}
       + s_\text{comp}V^\text{ldc}
       + s_\text{comp}^\text{off}V^\text{droop} \\
  0 &= -s_\text{frz} + \text{above}(V;\,V^\text{frz}) \\[0.5ex]
  0 &= -e_\text{RQ}
       + s_\text{ref}\left(V^\text{ref} - V^\text{meas}\right)
       + s_\text{ref}^\text{off}
         \left(k_\text{base}Q^\text{ref} - Q^\text{meas}\right) \\
  0 &= -e_\text{RQ}^\text{db}
       + \text{deadband2}
         \left(e_\text{RQ};\,D_\text{bd1},D_\text{bd2}\right) \\
  0 &= -e_\text{RQ}^\text{lim}
       + \text{clamp}
         \left(e_\text{RQ}^\text{db};\,e^{\min},e^{\max}\right) \\
  0 &= -Q^\text{PI}
       + \text{clamp}
         \left(K_\text{p}e_\text{RQ}^\text{lim}+x_Q^\text{PI};\,
               Q^{\min},Q^{\max}\right) \\
  0 &= -T_\text{fv}
         \left(k_\text{base}Q^\text{ext}-x_Q^\text{lag}\right)
       + T_\text{ft}
         \left(Q^\text{PI}-x_Q^\text{lag}\right) \\[0.5ex]
  0 &= -e_f
       + \text{deadband2}
         \left(f^\text{ref}-f;\,D_\text{bd1}^f,D_\text{bd2}^f\right) \\
  0 &= -e_P
       + k_\text{base}P_\text{plant}^\text{ref}
       - P^\text{meas}
       + D_\text{dn}\text{ramp}(e_f)
       - D_\text{up}\text{ramp}(-e_f) \\
  0 &= -e_P^\text{lim}
       + \text{clamp}\left(e_P;\,e_P^{\min},e_P^{\max}\right) \\
  0 &= -P^\text{PI}
       + \text{clamp}
         \left(K_\text{pg}e_P^\text{lim}+x_P^\text{PI};\,
               P^{\min},P^{\max}\right) \\
  0 &= -k_\text{base}P^\text{ext} + s_\text{freq}P^\text{ref}
\end{aligned}
```

CommonMath defines the [derived limiter functions](../../../../CommonMath.md#derived-functions)
and [ramp primitive](../../../../CommonMath.md#primitives) used above.

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_\text{r}, V_\text{i}
    &\leftarrow \text{regulated-bus voltage} \\
  I_\text{r}, I_\text{i}
    &\leftarrow \text{branch current} \\
  P^\text{br}, Q^\text{br}
    &\leftarrow \text{branch power} \\
  f
    &\leftarrow \text{frequency input} \\
  Q_0^\text{seed}
    &\leftarrow k_\text{base}Q_0^\text{ext} \\
  P_0^\text{seed}
    &\leftarrow k_\text{base}P_0^\text{ext}
\end{aligned}
```

### Internal Initialization

With initialization tolerance $\epsilon_0=10^{-10}$,
$\text{clamp}^{-1}(z;\ell,u)$ denotes a finite smooth-limiter input whose output
matches $z$ within $\epsilon_0$. It is available for $z\in[\ell,u]$,
including at both exact limits.

Subscript $0$ denotes initial values; all internal derivatives start at zero:

```math
\begin{aligned}
  V_0
    &= \sqrt{
       V_{\text{r},0}^2
       + V_{\text{i},0}^2
     } \\
  V_0^\text{ldc}
    &= \sqrt{
       \left(V_{\text{r},0}-R_c I_{\text{r},0}+X_c I_{\text{i},0}\right)^2
       + \left(V_{\text{i},0}-R_c I_{\text{i},0}-X_c I_{\text{r},0}\right)^2
     } \\
  V_0^\text{droop}
    &= V_0
       + K_c k_\text{base}Q_0^\text{br} \\
  V_0^\text{ctrl}
    &= s_\text{comp}V_0^\text{ldc}
       + s_\text{comp}^\text{off}V_0^\text{droop} \\
  V_0^\text{meas} &= V_0^\text{ctrl} \\
  Q_0^\text{meas} &= k_\text{base}Q_0^\text{br} \\
  P_0^\text{meas} &= k_\text{base}P_0^\text{br} \\
  s_{\text{frz},0}
    &= \text{above}\left(V_0;\,V^\text{frz}\right) \\
  e_{\text{RQ},0} &= 0 \\
  e_{\text{RQ},0}^\text{db}
    &= \text{deadband2}
       \left(e_{\text{RQ},0};\,D_\text{bd1},D_\text{bd2}\right) \\
  e_{\text{RQ},0}^\text{lim}
    &= \text{clamp}
       \left(e_{\text{RQ},0}^\text{db};\,e^{\min},e^{\max}\right) \\
  Q_0^\text{PI} &= Q_0^\text{seed} \\
  x_{Q,0}^\text{lag} &= Q_0^\text{PI} \\
  u_{Q,0}^\text{PI}
    &= \text{clamp}^{-1}
       \left(Q_0^\text{PI};\,Q^{\min},Q^{\max}\right) \\
  x_{Q,0}^\text{PI}
    &= u_{Q,0}^\text{PI}-K_\text{p}e_{\text{RQ},0}^\text{lim} \\
  e_{f,0}
    &= \text{deadband2}\left(0;\,D_\text{bd1}^f,D_\text{bd2}^f\right) \\
  P_0^\text{freq}
    &= D_\text{dn}\text{ramp}(e_{f,0})
       - D_\text{up}\text{ramp}(-e_{f,0}) \\
  e_{P,0} &= 0 \\
  e_{P,0}^\text{lim}
    &= \text{clamp}\left(e_{P,0};\,e_P^{\min},e_P^{\max}\right) \\
  P_0^\text{ref}
    &=
      \begin{cases}
        P_0^\text{seed} & s_\text{freq}=1 \\
        \text{clamp}
          \left(P_0^\text{meas};\,P^{\min},P^{\max}\right)
          & s_\text{freq}=0
      \end{cases} \\
  P_0^\text{PI} &= P_0^\text{ref} \\
  u_{P,0}^\text{PI}
    &=
      \begin{cases}
        \text{clamp}^{-1}
          \left(P_0^\text{PI};\,P^{\min},P^{\max}\right)
          & s_\text{freq}=1 \\
        P_0^\text{meas} & s_\text{freq}=0
      \end{cases} \\
  x_{P,0}^\text{PI}
    &= u_{P,0}^\text{PI}-K_\text{pg}e_{P,0}^\text{lim} \\
  P_0^\text{ext}
    &= \dfrac{s_\text{freq}}{k_\text{base}}P_0^\text{ref}
\end{aligned}
```

Initialization rejects an operating point when any of the following holds:

- a required bus or signal value, $Q_0^\text{seed}$, or a frequency-enabled
  $P_0^\text{seed}$ is not finite;
- $Q_0^\text{seed}\notin[Q^{\min},Q^{\max}]$, or
  $s_\text{freq}=1$ and
  $P_0^\text{seed}\notin[P^{\min},P^{\max}]$; or
- either PI antiwindup rate is nonfinite or has magnitude greater than
  $\epsilon_0$.

Every check resolves before any storage is written, so rejected initialization
leaves model and signal storage unchanged.

### Output Initialization

```math
\begin{aligned}
  f^\text{ref} &\leftarrow f_0 \\
  V^\text{ref} &\leftarrow V_0^\text{meas} \\
  Q^\text{ref} &\leftarrow Q_0^\text{br} \\
  P_\text{plant}^\text{ref}
    &\leftarrow \dfrac{P_0^\text{meas}-P_0^\text{freq}}{k_\text{base}}
\end{aligned}
```

## Monitorable Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`qext`          | [p.u.] | Reactive-power command output       | $Q^\text{ext}$ (system base)
`pext`          | [p.u.] | Active-power command output         | $P^\text{ext}$ (system base)
`vmeas`         | [p.u.] | Filtered regulated voltage          | $V^\text{meas}$
`qmeas`         | [p.u.] | Filtered reactive-power signal      | $Q^\text{meas}$ (component base)
`pmeas`         | [p.u.] | Filtered active-power signal        | $P^\text{meas}$ (component base)
