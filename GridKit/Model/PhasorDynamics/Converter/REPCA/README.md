# **Renewable Energy Plant Control Model (REPCA)**

REPCA is a WECC renewable energy plant control model for inverter-coupled
resources.

## Block Diagram

Standard REPCA block diagram.

![](../../../../../docs/Figures/PhasorDynamics/REPCA/diagram.png)

Figure 1: REPCA block diagram. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol                              | Units    | JSON        | Description                                             | Typical Value | Note
------------------------------------|----------|-------------|---------------------------------------------------------|---------------|------
$S^\mathrm{base}$                   | [MVA]    | `mva`       | REPCA component power base                              | 100.0         |
$s_{\mathrm{comp}}$                 | [binary] | `VcompFlag` | Voltage-compensation mode flag                          | 1.0           |
$s_{\mathrm{ref}}$                  | [binary] | `RefFlag`   | Reactive-loop reference flag                            | 1.0           |
$s_{\mathrm{freq}}$                 | [binary] | `Freqflag`  | Active-power control output flag                        | 0.0           |
$T_{\mathrm{fltr}}$                 | [sec]    | `Tfltr`     | Voltage and reactive-power measurement filter time constant | 0.05      |
$V^\mathrm{frz}$                    | [p.u.]   | `Vfrz`      | Regulated-voltage threshold below which the reactive PI state freezes | 0.7 |
$R_c$                               | [p.u.]   | `Rc`        | Line-drop compensation resistance                       | 0.0           |
$X_c$                               | [p.u.]   | `Xc`        | Line-drop compensation reactance                        | 0.0           |
$K_c$                               | [p.u.]   | `Kc`        | Reactive-current compensation gain                      | 1.0           |
$D_{\mathrm{bd1}}$                  | [p.u.]   | `dbdlow`    | Lower reactive-loop deadband threshold                  | 0.0           |
$D_{\mathrm{bd2}}$                  | [p.u.]   | `dbdupper`  | Upper reactive-loop deadband threshold                  | 0.0           |
$e^{\max}$                          | [p.u.]   | `emax`      | Maximum reactive-loop error limit                       | 1.0           |
$e^{\min}$                          | [p.u.]   | `emin`      | Minimum reactive-loop error limit                       | -1.0          |
$K_{\mathrm{p}}$                    | [p.u.]   | `Kp`        | Reactive-power controller proportional gain             | 10.0          |
$K_{\mathrm{i}}$                    | [p.u./s] | `Ki`        | Reactive-power controller integral gain                 | 10.0          |
$Q^{\max}$                          | [p.u.]   | `Qmax`      | Maximum reactive-power command                          | 1.0           |
$Q^{\min}$                          | [p.u.]   | `Qmin`      | Minimum reactive-power command                          | -1.0          |
$T_{\mathrm{ft}}$                   | [sec]    | `Tft`       | Reactive-command lead time constant                     | 0.0           |
$T_{\mathrm{fv}}$                   | [sec]    | `Tfv`       | Reactive-command lag time constant                      | 3.0           |
$T_{\mathrm{p}}$                    | [sec]    | `Tp`        | Active-power measurement filter time constant           | 0.0           |
$D_{\mathrm{bd1}}^f$                | [p.u.]   | `fdbd1`     | Lower frequency-error deadband threshold                | 0.0           |
$D_{\mathrm{bd2}}^f$                | [p.u.]   | `fdbd2`     | Upper frequency-error deadband threshold                | 0.0           |
$D_{\mathrm{dn}}$                   | [p.u.]   | `Ddn`       | Down-frequency droop response gain                      | 20.0          |
$D_{\mathrm{up}}$                   | [p.u.]   | `Dup`       | Up-frequency droop response gain                        | 0.0           |
$e_P^{\max}$                        | [p.u.]   | `femax`     | Maximum active-power error limit                        | 1.0           |
$e_P^{\min}$                        | [p.u.]   | `femin`     | Minimum active-power error limit                        | -1.0          |
$K_{\mathrm{pg}}$                   | [p.u.]   | `Kpg`       | Active-power controller proportional gain               | 10.0          |
$K_{\mathrm{ig}}$                   | [p.u./s] | `Kig`       | Active-power controller integral gain                   | 10.0          |
$P^{\max}$                          | [p.u.]   | `Pmax`      | Maximum active-power command                            | 2.0           |
$P^{\min}$                          | [p.u.]   | `Pmin`      | Minimum active-power command                            | 0.0           |
$T_{\mathrm{lag}}$                  | [sec]    | `Tlag`      | Active-power command lag time constant                  | 3.0           |

### Parameter Validation

Invalid REPCA parameter sets are rejected by the following checks. Let $\epsilon_T=10^{-3}$.

```math
\begin{aligned}
  T &\leftarrow \max\!\left(T, \epsilon_T\right)
    \quad T\in\{T_{\mathrm{fltr}},T_{\mathrm{p}},T_{\mathrm{lag}}\} \\
  S^\mathrm{base}
    &> 0 \\
  s_{\mathrm{comp}}, s_{\mathrm{ref}}, s_{\mathrm{freq}}
    &\in \{0,1\} \\
  T_{\mathrm{fv}}
    &> 0 \\
  D_{\mathrm{bd1}}
    &\le 0 \le D_{\mathrm{bd2}} \\
  e^{\min}
    &\le 0 \le e^{\max} \\
  Q^{\min}
    &\le Q^{\max} \\
  D_{\mathrm{bd1}}^f
    &\le 0 \le D_{\mathrm{bd2}}^f \\
  e_P^{\min}
    &\le 0 \le e_P^{\max} \\
  P^{\min}
    &\le P^{\max}
\end{aligned}
```

### Model Derived Parameters

```math
\begin{aligned}
  s_{\mathrm{comp}}^\mathrm{off} &= 1 - s_{\mathrm{comp}} \\
  s_{\mathrm{ref}}^\mathrm{off} &= 1 - s_{\mathrm{ref}} \\
  k_{\mathrm{base}} &= \dfrac{S^\mathrm{sys}}{S^\mathrm{base}}
\end{aligned}
```

## Model Ports

Name        | Port   | Init    | Description
------------|--------|---------|------
`bus`       | Bus    | Known   | Regulated bus voltage
`ibranchr`  | Input  | Known   | Branch current real component
`ibranchi`  | Input  | Known   | Branch current imaginary component
`pbranch`   | Input  | Known   | Branch active power
`qbranch`   | Input  | Known   | Branch reactive power
`freq`      | Input  | Known   | Frequency input
`freqref`   | Input  | Unknown | Frequency reference
`vref`      | Input  | Unknown | Voltage-control reference
`qref`      | Input  | Unknown | Reactive-power reference
`pplantref` | Input  | Unknown | Plant active-power reference
`qext`      | Output | Known   | Reactive-power command output
`pext`      | Output | Known   | Active-power command output

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                         | Note
------------------------|--------|-------------------------------------|------
$V^\mathrm{meas}$       | [p.u.] | Filtered regulated voltage          | State 1 in Fig. 1; source label: `Vmeas`
$Q^\mathrm{meas}$       | [p.u.] | Filtered reactive-power signal      | State 2 in Fig. 1; source label: `Qmeas`
$x_Q^\mathrm{PI}$       | [p.u.] | Reactive PI controller state        | State 3 in Fig. 1; source label: `Reactive PI`
$x_Q^\mathrm{lag}$      | [p.u.] | Reactive-command lead-lag state     | State 4 in Fig. 1; source label: `Qext`
$P^\mathrm{meas}$       | [p.u.] | Filtered active-power signal        | State 5 in Fig. 1; source label: `Pmeas`
$x_P^\mathrm{PI}$       | [p.u.] | Active-power PI controller state    | State 6 in Fig. 1; source label: `Power PI`
$P^\mathrm{ref}$        | [p.u.] | Active-power command lag state      | State 7 in Fig. 1; source label: `Pref`

#### Algebraic

Symbol                          | Units  | Description                         | Note
--------------------------------|--------|-------------------------------------|------
$V$                             | [p.u.] | Regulated-bus voltage magnitude     |
$V^\mathrm{ldc}$                  | [p.u.] | Line-drop compensated voltage magnitude |
$V^\mathrm{droop}$                | [p.u.] | Reactive-droop-compensated voltage |
$V^\mathrm{ctrl}$                 | [p.u.] | Selected voltage-measurement input  |
$s_{\mathrm{frz}}$              | [binary] | Reactive-PI voltage-enable indicator |
$e_{\mathrm{RQ}}$               | [p.u.] | Selected reactive-loop error        |
$e_{\mathrm{RQ}}^\mathrm{db}$     | [p.u.] | Deadbanded reactive-loop error      |
$e_{\mathrm{RQ}}^\mathrm{lim}$    | [p.u.] | Limited reactive-loop error         |
$Q^\mathrm{PI}$                  | [p.u.] | Reactive PI output                  |
$Q^\mathrm{ext}$                | [p.u.] | Reactive-power command output       | System base
$e_f$                           | [p.u.] | Frequency error after deadband      |
$e_P$                           | [p.u.] | Active-power control error          |
$e_P^\mathrm{lim}$                | [p.u.] | Limited active-power control error  |
$P^\mathrm{PI}$                  | [p.u.] | Active-power PI output              |
$P^\mathrm{ext}$                | [p.u.] | Active-power command output         | System base

### External Variables

#### Differential
None.

#### Algebraic

Symbol                               | Units  | Type   | Description                       | Note
-------------------------------------|--------|--------|-----------------------------------|------
$V_{\mathrm{r}}$                     | [p.u.] | Known   | Regulated-bus voltage, real component | Source label: `Vreg`
$V_{\mathrm{i}}$                     | [p.u.] | Known   | Regulated-bus voltage, imaginary component | Source label: `Vreg`
$I_{\mathrm{r}}$                     | [p.u.] | Known   | Branch current real component       | Source label: `Ibranch`
$I_{\mathrm{i}}$                     | [p.u.] | Known   | Branch current imaginary component  | Source label: `Ibranch`
$P^\mathrm{br}$                       | [p.u.] | Known   | Branch active power                 |
$Q^\mathrm{br}$                       | [p.u.] | Known   | Branch reactive power               |
$f$                                  | [p.u.] | Known   | Frequency input                    | Source label: `Freq`
$f^\mathrm{ref}$                      | [p.u.] | Unknown | Frequency reference                | Optional signal port `freqref`; source label: `Freq_ref`
$V^\mathrm{ref}$                      | [p.u.] | Unknown | Voltage-control reference          | Optional signal port `vref`
$Q^\mathrm{ref}$                      | [p.u.] | Unknown | Reactive-power reference           | Optional signal port `qref`; system base
$P_\mathrm{plant}^\mathrm{ref}$     | [p.u.] | Unknown | Plant active-power reference       | Optional signal port `pplantref`; system base

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{V}^\mathrm{meas}
    + \dfrac{1}{T_{\mathrm{fltr}}}
      \left(V^\mathrm{ctrl} - V^\mathrm{meas}\right) \\
  0 &=
    -\dot{Q}^\mathrm{meas}
    + \dfrac{1}{T_{\mathrm{fltr}}}
      \left(k_{\mathrm{base}}Q^\mathrm{br} - Q^\mathrm{meas}\right) \\
  0 &=
    -\dot{x}_Q^\mathrm{PI}
    + s_{\mathrm{frz}}
    \text{antiwindup}\!\left(
      Q^\mathrm{PI},
      K_{\mathrm{i}}e_{\mathrm{RQ}}^\mathrm{lim},
      Q^{\min},
      Q^{\max}
    \right) \\
  0 &=
    -\dot{x}_Q^\mathrm{lag}
    + \dfrac{1}{T_{\mathrm{fv}}}
      \left(Q^\mathrm{PI} - x_Q^\mathrm{lag}\right) \\
  0 &=
    -\dot{P}^\mathrm{meas}
    + \dfrac{1}{T_{\mathrm{p}}}
      \left(k_{\mathrm{base}}P^\mathrm{br} - P^\mathrm{meas}\right) \\
  0 &=
    -\dot{x}_P^\mathrm{PI}
    + \text{antiwindup}\!\left(
      P^\mathrm{PI},
      K_{\mathrm{ig}}e_P^\mathrm{lim},
      P^{\min},
      P^{\max}
    \right) \\
  0 &=
    -\dot{P}^\mathrm{ref}
    + \dfrac{1}{T_{\mathrm{lag}}}
      \left(P^\mathrm{PI} - P^\mathrm{ref}\right)
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#anti-windup-indicator) target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -V^2
       + V_{\mathrm{r}}^2
       + V_{\mathrm{i}}^2 \\
  0 &=
    -\left(V^\mathrm{ldc}\right)^2
    + \left(V_{\mathrm{r}} - R_c I_{\mathrm{r}} + X_c I_{\mathrm{i}}\right)^2
    + \left(V_{\mathrm{i}} - R_c I_{\mathrm{i}} - X_c I_{\mathrm{r}}\right)^2 \\
  0 &= -V^\mathrm{droop} + V + K_c k_{\mathrm{base}}Q^\mathrm{br} \\
  0 &= -V^\mathrm{ctrl} + s_{\mathrm{comp}}V^\mathrm{ldc} + s_{\mathrm{comp}}^\mathrm{off}V^\mathrm{droop} \\
  0 &= -s_{\mathrm{frz}} + \text{above}(V, V^\mathrm{frz}) \\[0.5ex]
  0 &= -e_{\mathrm{RQ}}
       + s_{\mathrm{ref}}\left(V^\mathrm{ref} - V^\mathrm{meas}\right)
       + s_{\mathrm{ref}}^\mathrm{off}\left(k_{\mathrm{base}}Q^\mathrm{ref} - Q^\mathrm{meas}\right) \\
  0 &= -e_{\mathrm{RQ}}^\mathrm{db}
       + \text{deadband2}(e_{\mathrm{RQ}}, D_{\mathrm{bd1}}, D_{\mathrm{bd2}}) \\
  0 &= -e_{\mathrm{RQ}}^\mathrm{lim}
       + \text{clamp}(e_{\mathrm{RQ}}^\mathrm{db}, e^{\min}, e^{\max}) \\
  0 &= -Q^\mathrm{PI} + \text{clamp}(K_{\mathrm{p}} e_{\mathrm{RQ}}^\mathrm{lim} + x_Q^\mathrm{PI}, Q^{\min}, Q^{\max}) \\
  0 &= -T_{\mathrm{fv}}(k_{\mathrm{base}}Q^\mathrm{ext} - x_Q^\mathrm{lag})
       + T_{\mathrm{ft}}(Q^\mathrm{PI} - x_Q^\mathrm{lag}) \\[0.5ex]
  0 &= -e_f + \text{deadband2}(f^\mathrm{ref} - f, D_{\mathrm{bd1}}^f, D_{\mathrm{bd2}}^f) \\
  0 &= -e_P
       + k_{\mathrm{base}}P_\mathrm{plant}^\mathrm{ref}
       - P^\mathrm{meas}
       + D_{\mathrm{dn}}\rho(e_f)
       - D_{\mathrm{up}}\rho(-e_f) \\
  0 &= -e_P^\mathrm{lim} + \text{clamp}(e_P, e_P^{\min}, e_P^{\max}) \\
  0 &= -P^\mathrm{PI} + \text{clamp}(K_{\mathrm{pg}} e_P^\mathrm{lim} + x_P^\mathrm{PI}, P^{\min}, P^{\max}) \\
  0 &= -k_{\mathrm{base}}P^\mathrm{ext} + s_{\mathrm{freq}}P^\mathrm{ref}
\end{aligned}
```

CommonMath defines [above](../../../../CommonMath.md#derived-functions),
[clamp](../../../../CommonMath.md#derived-functions),
[deadband2](../../../../CommonMath.md#derived-functions), and
[ramp](../../../../CommonMath.md#primitives) $\rho$.

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_{\mathrm{r}}, V_{\mathrm{i}}
    &\leftarrow \text{regulated-bus voltage} \\
  I_{\mathrm{r}}, I_{\mathrm{i}}
    &\leftarrow \text{branch current} \\
  P^\mathrm{br}, Q^\mathrm{br}
    &\leftarrow \text{branch power} \\
  f
    &\leftarrow \text{frequency input} \\
  Q^\mathrm{ext}
    &\leftarrow \text{reactive-power command start} \\
  P^\mathrm{ext}
    &\leftarrow \text{active-power command start}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  V_0
    &= \sqrt{
       V_{\mathrm{r},0}^2
       + V_{\mathrm{i},0}^2
     } \\
  V_0^\mathrm{ldc}
    &= \sqrt{
       \left(V_{\mathrm{r},0} - R_c I_{\mathrm{r},0} + X_c I_{\mathrm{i},0}\right)^2
       + \left(V_{\mathrm{i},0} - R_c I_{\mathrm{i},0} - X_c I_{\mathrm{r},0}\right)^2
     } \\
  V_0^\mathrm{droop}
    &= V_0
       + K_c k_{\mathrm{base}}Q_0^\mathrm{br} \\
  V_0^\mathrm{ctrl}
    &= s_{\mathrm{comp}}V_0^\mathrm{ldc}
       + s_{\mathrm{comp}}^\mathrm{off}V_0^\mathrm{droop} \\
  e_{f,0} &= 0 \\
  P_0^\mathrm{freq}
    &= D_{\mathrm{dn}}\rho(e_{f,0})
       - D_{\mathrm{up}}\rho(-e_{f,0}) \\
  V_0^\mathrm{meas} &= V_0^\mathrm{ctrl} \\
  Q_0^\mathrm{meas} &= k_{\mathrm{base}}Q_0^\mathrm{br} \\
  P_0^\mathrm{meas} &= k_{\mathrm{base}}P_0^\mathrm{br} \\
  s_{\mathrm{frz},0}
    &= \text{above}\left(V_0, V^\mathrm{frz}\right) \\
  e_{\mathrm{RQ},0} &= 0 \\
  e_{\mathrm{RQ},0}^\mathrm{db}
    &= \text{deadband2}\left(e_{\mathrm{RQ},0}, D_{\mathrm{bd1}}, D_{\mathrm{bd2}}\right) \\
  e_{\mathrm{RQ},0}^\mathrm{lim}
    &= \text{clamp}\left(e_{\mathrm{RQ},0}^\mathrm{db}, e^{\min}, e^{\max}\right) \\
  Q_0^\mathrm{PI} &= k_{\mathrm{base}}Q_0^\mathrm{ext} \\
  x_{Q,0}^\mathrm{lag} &= k_{\mathrm{base}}Q_0^\mathrm{ext} \\
  x_{Q,0}^\mathrm{PI}
    &= k_{\mathrm{base}}Q_0^\mathrm{ext} - K_{\mathrm{p}}e_{\mathrm{RQ},0}^\mathrm{lim} \\
  e_{P,0} &= 0 \\
  e_{P,0}^\mathrm{lim}
    &= \text{clamp}\left(e_{P,0}, e_P^{\min}, e_P^{\max}\right) \\
  P_0^\mathrm{ref}
    &=
      \begin{cases}
        k_{\mathrm{base}}P_0^\mathrm{ext} & s_{\mathrm{freq}} = 1 \\
        \text{clamp}
          \left(k_{\mathrm{base}}P_0^\mathrm{br}, P^{\min}, P^{\max}\right)
          & s_{\mathrm{freq}} = 0
      \end{cases} \\
  P_0^\mathrm{PI} &= P_0^\mathrm{ref} \\
  x_{P,0}^\mathrm{PI}
    &= P_0^\mathrm{ref} - K_{\mathrm{pg}}e_{P,0}^\mathrm{lim}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
  f^\mathrm{ref} &\leftarrow f_0 \\
  V^\mathrm{ref} &\leftarrow V_0^\mathrm{meas} \\
  Q^\mathrm{ref} &\leftarrow \dfrac{Q_0^\mathrm{meas}}{k_{\mathrm{base}}} \\
  P_{\mathrm{plant}}^\mathrm{ref}
    &\leftarrow \dfrac{P_0^\mathrm{meas} - P_0^\mathrm{freq}}{k_{\mathrm{base}}}
\end{aligned}
```

Initialization rejects nonzero reactive-power or active-power PI antiwindup
rates.

REPCA writes the resolved references to attached `freqref`, `vref`, `qref`,
and `pplantref` signal inputs. If no controller is connected, those values are
used as constant reference inputs.

## Monitorable Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`qext`          | [p.u.] | Reactive-power command output       | $Q^\mathrm{ext}$ (system base)
`pext`          | [p.u.] | Active-power command output         | $P^\mathrm{ext}$ (system base)
`vmeas`         | [p.u.] | Filtered regulated voltage          | $V^\mathrm{meas}$
`qmeas`         | [p.u.] | Filtered reactive-power signal      | $Q^\mathrm{meas}$ (component base)
`pmeas`         | [p.u.] | Filtered active-power signal        | $P^\mathrm{meas}$ (component base)
