# **Renewable Energy Plant Control Model (REPCA)**

REPCA is a WECC renewable energy plant control model for inverter-coupled
resources. In GridKit it is represented as a plant-level signal-control model
that computes active- and reactive-power commands for downstream electrical
control models.

## Block Diagram

Standard REPCA block diagram.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_REPCA_Diagram.png">

  Figure 1: REPCA block diagram. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                              | Units    | JSON        | Description                                             | Typical Value | Note
------------------------------------|----------|-------------|---------------------------------------------------------|---------------|------
$S^{\mathrm{base}}$                 | [MVA]    | `mva`       | REPCA component power base                              | 0.0           | Block name: `MVABase`; zero uses system base
$s_{\mathrm{comp}}$                 | [binary] | `VcompFlag` | Voltage-compensation mode flag                          | 1.0           | 1 = line-drop compensation, 0 = reactive-droop compensation
$s_{\mathrm{ref}}$                  | [binary] | `RefFlag`   | Reactive-loop reference flag                            | 1.0           | 1 = voltage control, 0 = reactive-power control
$s_{\mathrm{freq}}$                 | [binary] | `Freqflag`  | Active-power control output flag                        | 0.0           | 1 = use active-power loop output, 0 = output zero
$T_{\mathrm{fltr}}$                 | [sec]    | `Tfltr`     | Voltage and reactive-power measurement filter time constant | 0.05      | If zero, $V_{\mathrm{meas}}$ and $Q_{\mathrm{meas}}$ are algebraic
$T_{\mathrm{ft}}$                   | [sec]    | `Tft`       | Reactive-command lead time constant                     | 0.0           |
$T_{\mathrm{fv}}$                   | [sec]    | `Tfv`       | Reactive-command lag time constant                      | 3.0           |
$T_{\mathrm{p}}$                    | [sec]    | `Tp`        | Active-power measurement filter time constant           | 0.0           | If zero, $P_{\mathrm{meas}}$ is algebraic
$T_{\mathrm{lag}}$                  | [sec]    | `Tlag`      | Active-power command lag time constant                  | 3.0           | If zero, $P_{\mathrm{ref}}$ is algebraic
$V_{\mathrm{frz}}$                  | [p.u.]   | `Vfrz`      | Regulated-voltage threshold below which the reactive PI state freezes | 0.7 |
$R_c$                               | [p.u.]   | `Rc`        | Line-drop compensation resistance                       | 0.0           |
$X_c$                               | [p.u.]   | `Xc`        | Line-drop compensation reactance                        | 0.0           |
$K_c$                               | [p.u.]   | `Kc`        | Reactive-current compensation gain                      | 1.0           |
$D_{\mathrm{bd1}}$                  | [p.u.]   | `dbdlow`    | Lower reactive-loop deadband threshold                  | 0.0           | REPCA1/REPCTA1 key; negative of `dbd` for symmetric REPC_A data
$D_{\mathrm{bd2}}$                  | [p.u.]   | `dbdupper`  | Upper reactive-loop deadband threshold                  | 0.0           | REPCA1/REPCTA1 key; `dbd` for symmetric REPC_A data
$e^{\max}$                          | [p.u.]   | `emax`      | Maximum reactive-loop error limit                       | 1.0           |
$e^{\min}$                          | [p.u.]   | `emin`      | Minimum reactive-loop error limit                       | -1.0          |
$K_{\mathrm{p}}$                    | [p.u.]   | `Kp`        | Reactive-power controller proportional gain             | 10.0          |
$K_{\mathrm{i}}$                    | [p.u./s] | `Ki`        | Reactive-power controller integral gain                 | 10.0          |
$Q^{\max}$                          | [p.u.]   | `Qmax`      | Maximum reactive-power command                          | 1.0           |
$Q^{\min}$                          | [p.u.]   | `Qmin`      | Minimum reactive-power command                          | -1.0          |
$D_{f,\mathrm{bd1}}$                | [p.u.]   | `fdbd1`     | Lower frequency-error deadband threshold                | 0.0           |
$D_{f,\mathrm{bd2}}$                | [p.u.]   | `fdbd2`     | Upper frequency-error deadband threshold                | 0.0           |
$D_{\mathrm{dn}}$                   | [p.u.]   | `Ddn`       | Down-frequency droop response gain                      | 20.0          |
$D_{\mathrm{up}}$                   | [p.u.]   | `Dup`       | Up-frequency droop response gain                        | 0.0           |
$e_P^{\max}$                        | [p.u.]   | `femax`     | Maximum active-power error limit                        | 1.0           |
$e_P^{\min}$                        | [p.u.]   | `femin`     | Minimum active-power error limit                        | -1.0          |
$K_{\mathrm{pg}}$                   | [p.u.]   | `Kpg`       | Active-power controller proportional gain               | 10.0          |
$K_{\mathrm{ig}}$                   | [p.u./s] | `Kig`       | Active-power controller integral gain                   | 10.0          |
$P^{\max}$                          | [p.u.]   | `Pmax`      | Maximum active-power command                            | 2.0           |
$P^{\min}$                          | [p.u.]   | `Pmin`      | Minimum active-power command                            | 0.0           |

### Parameter Validation

Invalid REPCA parameter sets are rejected by the following checks.

The required checks are:

```math
\begin{aligned}
  &S^{\mathrm{base}} \ge 0 \\
  &s_{\mathrm{comp}}, s_{\mathrm{ref}}, s_{\mathrm{freq}} \in \{0,1\} \\
  &T_{\mathrm{fltr}}, T_{\mathrm{p}}, T_{\mathrm{lag}} \ge 0 \\
  &T_{\mathrm{ft}} \ge 0,\quad T_{\mathrm{fv}} > 0 \\
  &V_{\mathrm{frz}} \ge 0 \\
  &D_{\mathrm{bd1}} \le 0 \le D_{\mathrm{bd2}} \\
  &e^{\min} \le 0 \le e^{\max} \\
  &Q^{\min} \le Q^{\max} \\
  &D_{f,\mathrm{bd1}} \le 0 \le D_{f,\mathrm{bd2}} \\
  &D_{\mathrm{dn}}, D_{\mathrm{up}} \ge 0 \\
  &e_P^{\min} \le 0 \le e_P^{\max} \\
  &P^{\min} \le P^{\max}
\end{aligned}
```

### Model Derived Parameters

The off-mode flag complements are:

```math
\begin{aligned}
  s_{\mathrm{comp}}^{\mathrm{off}} &= 1 - s_{\mathrm{comp}} \\
  s_{\mathrm{ref}}^{\mathrm{off}} &= 1 - s_{\mathrm{ref}}
\end{aligned}
```

The branch-power base conversion factor is:

```math
\begin{aligned}
  k_{\mathrm{base}} &=
  \begin{cases}
    1, & S^{\mathrm{base}} = 0 \\
    \dfrac{S^{\mathrm{sys}}}{S^{\mathrm{base}}}, & S^{\mathrm{base}} > 0
  \end{cases}
\end{aligned}
```

Here $S^{\mathrm{sys}}$ is the network system power base. Branch power inputs
are on system base and are converted to component base with
$k_{\mathrm{base}}$. Branch current inputs are already on component base.

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                         | Note
------------------------|--------|-------------------------------------|------
$V_{\mathrm{meas}}$     | [p.u.] | Filtered regulated voltage          | State 1 in Fig. 1; source label: `Vmeas`; algebraic when $T_{\mathrm{fltr}} = 0$
$Q_{\mathrm{meas}}$     | [p.u.] | Filtered reactive-power signal      | State 2 in Fig. 1; source label: `Qmeas`; algebraic when $T_{\mathrm{fltr}} = 0$
$x_Q$                   | [p.u.] | Reactive PI controller state        | State 3 in Fig. 1; source label: `Reactive PI`
$x_{\mathrm{Qext}}$     | [p.u.] | Reactive-command lead-lag state     | State 4 in Fig. 1; source label: `Qext`
$P_{\mathrm{meas}}$     | [p.u.] | Filtered active-power signal        | State 5 in Fig. 1; source label: `Pmeas`; algebraic when $T_{\mathrm{p}} = 0$
$x_P$                   | [p.u.] | Active-power PI controller state    | State 6 in Fig. 1; source label: `Power PI`
$P_{\mathrm{ref}}$      | [p.u.] | Active-power command lag state      | State 7 in Fig. 1; source label: `Pref`; algebraic when $T_{\mathrm{lag}} = 0$

#### Algebraic

Symbol                          | Units  | Description                         | Note
--------------------------------|--------|-------------------------------------|------
$V_{\mathrm{reg}}$              | [p.u.] | Regulated-bus voltage magnitude     | From regulated-bus phasor
$V_{\mathrm{ldc}}$              | [p.u.] | Line-drop compensated voltage magnitude | Selected when $s_{\mathrm{comp}}=1$
$V_{\mathrm{droop}}$            | [p.u.] | Reactive-droop-compensated voltage | Selected when $s_{\mathrm{comp}}=0$
$V_{\mathrm{ctrl}}$             | [p.u.] | Selected voltage-measurement input  | Input to $V_{\mathrm{meas}}$ filter
$s_{\mathrm{frz}}$              | [binary] | Reactive-PI voltage-enable indicator | 1 when $V_{\mathrm{reg}} > V_{\mathrm{frz}}$
$e_{\mathrm{RQ}}$               | [p.u.] | Selected reactive-loop error        | Chosen by $s_{\mathrm{ref}}$
$e_{\mathrm{RQ}}^{\mathrm{db}}$  | [p.u.] | Deadbanded reactive-loop error      | Defined by CommonMath `deadband2`
$e_{\mathrm{RQ}}^{\mathrm{lim}}$ | [p.u.] | Limited reactive-loop error         | Feeds reactive PI
$Q_{\mathrm{PI}}$               | [p.u.] | Reactive PI output                  | Limited by $Q^{\min}$ and $Q^{\max}$
$Q_{\mathrm{ext}}$              | [p.u.] | Reactive-power command output       | Sent to downstream electrical control
$e_f$                           | [p.u.] | Frequency error after deadband      | From $f_{\mathrm{ref}} - f$
$e_P$                           | [p.u.] | Active-power control error          | Plant reference minus $P_{\mathrm{meas}}$ plus frequency droop correction
$e_P^{\mathrm{lim}}$            | [p.u.] | Limited active-power control error  | Feeds active-power PI
$P_{\mathrm{PI}}$               | [p.u.] | Active-power PI output              | Limited by $P^{\min}$ and $P^{\max}$
$P_{\mathrm{ext}}$              | [p.u.] | Active-power command output         | Sent to downstream electrical control when $s_{\mathrm{freq}}=1$

### External Variables

#### Differential
None.

#### Algebraic

Symbol                               | Units  | Description                       | Note
-------------------------------------|--------|-----------------------------------|------
$V_{\mathrm{reg,r}}$                 | [p.u.] | Regulated-bus voltage, real component | Source label: `Vreg`
$V_{\mathrm{reg,i}}$                 | [p.u.] | Regulated-bus voltage, imaginary component | Source label: `Vreg`
$I_{\mathrm{br,r}}$                  | [p.u.] | Branch current real component       | Component base; source label: `Ibranch`; required
$I_{\mathrm{br,i}}$                  | [p.u.] | Branch current imaginary component  | Component base; source label: `Ibranch`; required
$P_{\mathrm{br}}$                    | [p.u.] | Branch active power                 | System base; required
$Q_{\mathrm{br}}$                    | [p.u.] | Branch reactive power               | System base; required
$V_{\mathrm{ref}}$                   | [p.u.] | Voltage-control reference          | Optional, defaults to initialized constant
$Q_{\mathrm{ref}}$                   | [p.u.] | Reactive-power reference           | Component base; optional, defaults to initialized constant
$P_{\mathrm{plant}}^{\mathrm{ref}}$  | [p.u.] | Plant active-power reference       | Component base; optional, defaults to initialized constant
$f$                                  | [p.u.] | Frequency input                    | Source label: `Freq`; optional, defaults to zero
$f_{\mathrm{ref}}$                   | [p.u.] | Frequency reference                | Source label: `Freq_ref`; optional, defaults to zero

## Model Equations

### Differential Equations

The measurement filters and active-power output lag are written in descriptor form; when $T_{\mathrm{fltr}} = 0$, $T_{\mathrm{p}} = 0$, or $T_{\mathrm{lag}} = 0$, the corresponding residual becomes algebraic. The reactive-command lead-lag denominator requires $T_{\mathrm{fv}} > 0$.

```math
\begin{aligned}
  0 &= -T_{\mathrm{fltr}}\dot V_{\mathrm{meas}} - V_{\mathrm{meas}} + V_{\mathrm{ctrl}} \\
  0 &= -T_{\mathrm{fltr}}\dot Q_{\mathrm{meas}} - Q_{\mathrm{meas}} + k_{\mathrm{base}}Q_{\mathrm{br}} \\
  0 &=
    -\dot x_Q
    + s_{\mathrm{frz}}
    \text{antiwindup}\!\left(
      Q_{\mathrm{PI}},
      K_{\mathrm{i}}e_{\mathrm{RQ}}^{\mathrm{lim}},
      Q^{\min},
      Q^{\max}
    \right) \\
  0 &= -T_{\mathrm{fv}}\dot x_{\mathrm{Qext}} - x_{\mathrm{Qext}} + Q_{\mathrm{PI}} \\
  0 &= -T_{\mathrm{p}}\dot P_{\mathrm{meas}} - P_{\mathrm{meas}} + k_{\mathrm{base}}P_{\mathrm{br}} \\
  0 &=
    -\dot x_P
    + \text{antiwindup}\!\left(
      P_{\mathrm{PI}},
      K_{\mathrm{ig}}e_P^{\mathrm{lim}},
      P^{\min},
      P^{\max}
    \right) \\
  0 &= -T_{\mathrm{lag}}\dot P_{\mathrm{ref}} - P_{\mathrm{ref}} + P_{\mathrm{PI}}
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#anti-windup-indicator) target and smooth approximation.

### Algebraic Equations

The algebraic targets use CommonMath helper notation where applicable:

```math
\begin{aligned}
  0 &= -V_{\mathrm{reg}}^2 + V_{\mathrm{reg,r}}^2 + V_{\mathrm{reg,i}}^2 \\
  0 &=
    -V_{\mathrm{ldc}}^2
    + (V_{\mathrm{reg,r}} - R_c I_{\mathrm{br,r}} + X_c I_{\mathrm{br,i}})^2
    + (V_{\mathrm{reg,i}} - R_c I_{\mathrm{br,i}} - X_c I_{\mathrm{br,r}})^2 \\
  0 &= -V_{\mathrm{droop}} + V_{\mathrm{reg}} + K_c k_{\mathrm{base}}Q_{\mathrm{br}} \\
  0 &= -V_{\mathrm{ctrl}} + s_{\mathrm{comp}}V_{\mathrm{ldc}} + s_{\mathrm{comp}}^{\mathrm{off}}V_{\mathrm{droop}} \\
  0 &= -s_{\mathrm{frz}} + \text{above}(V_{\mathrm{reg}}, V_{\mathrm{frz}}) \\[0.5ex]
  0 &= -e_{\mathrm{RQ}}
       + s_{\mathrm{ref}}\left(V_{\mathrm{ref}} - V_{\mathrm{meas}}\right)
       + s_{\mathrm{ref}}^{\mathrm{off}}\left(Q_{\mathrm{ref}} - Q_{\mathrm{meas}}\right) \\
  0 &= -e_{\mathrm{RQ}}^{\mathrm{db}}
       + \text{deadband2}(e_{\mathrm{RQ}}, D_{\mathrm{bd1}}, D_{\mathrm{bd2}}) \\
  0 &= -e_{\mathrm{RQ}}^{\mathrm{lim}}
       + \text{clamp}(e_{\mathrm{RQ}}^{\mathrm{db}}, e^{\min}, e^{\max}) \\
  0 &= -Q_{\mathrm{PI}} + \text{clamp}(K_{\mathrm{p}} e_{\mathrm{RQ}}^{\mathrm{lim}} + x_Q, Q^{\min}, Q^{\max}) \\
  0 &= -T_{\mathrm{fv}}(Q_{\mathrm{ext}} - x_{\mathrm{Qext}})
       + T_{\mathrm{ft}}(Q_{\mathrm{PI}} - x_{\mathrm{Qext}}) \\[0.5ex]
  0 &= -e_f + \text{deadband2}(f_{\mathrm{ref}} - f, D_{f,\mathrm{bd1}}, D_{f,\mathrm{bd2}}) \\
  0 &= -e_P
       + P_{\mathrm{plant}}^{\mathrm{ref}}
       - P_{\mathrm{meas}}
       + D_{\mathrm{dn}}\rho(e_f)
       - D_{\mathrm{up}}\rho(-e_f) \\
  0 &= -e_P^{\mathrm{lim}} + \text{clamp}(e_P, e_P^{\min}, e_P^{\max}) \\
  0 &= -P_{\mathrm{PI}} + \text{clamp}(K_{\mathrm{pg}} e_P^{\mathrm{lim}} + x_P, P^{\min}, P^{\max}) \\
  0 &= -P_{\mathrm{ext}} + s_{\mathrm{freq}}P_{\mathrm{ref}}
\end{aligned}
```

The $V_{\mathrm{reg}}$ and $V_{\mathrm{ldc}}$ variables use nonnegative branches of squared algebraic residuals. The exact frequency-droop target is $D_{\mathrm{dn}}e_f$ for positive deadbanded frequency error, $D_{\mathrm{up}}e_f$ for negative deadbanded frequency error, and zero in the deadband. The $\rho$ form above is the smooth CommonMath representation of that split.

CommonMath defines the helper targets and smooth approximations for [above, clamp, and deadband2](../../../../CommonMath.md#derived-functions). The frequency split uses the primitive [ramp](../../../../CommonMath.md#primitives) $\rho$.

## Initialization

Initialization is performed by evaluating the steady-state residuals in dependency order. Let subscript $0$ denote initial values and set all internal derivatives to zero. REPCA reads branch power and branch current from connected signal ports:

```math
\begin{aligned}
  P_{\mathrm{br},0} &= pbranch_0 \\
  Q_{\mathrm{br},0} &= qbranch_0 \\
  I_{\mathrm{br,r},0} &= ibranchr_0 \\
  I_{\mathrm{br,i},0} &= ibranchi_0
\end{aligned}
```

Then convert branch power to component base and evaluate the voltage-compensation signals:

```math
\begin{aligned}
  V_{\mathrm{reg},0} &= \sqrt{V_{\mathrm{reg,r},0}^2 + V_{\mathrm{reg,i},0}^2} \\
  V_{\mathrm{ldc},0}
    &= \sqrt{
       (V_{\mathrm{reg,r},0} - R_c I_{\mathrm{br,r},0} + X_c I_{\mathrm{br,i},0})^2
       + (V_{\mathrm{reg,i},0} - R_c I_{\mathrm{br,i},0} - X_c I_{\mathrm{br,r},0})^2
     } \\
  V_{\mathrm{droop},0} &= V_{\mathrm{reg},0} + K_c k_{\mathrm{base}}Q_{\mathrm{br},0} \\
  V_{\mathrm{ctrl},0} &= s_{\mathrm{comp}}V_{\mathrm{ldc},0} + s_{\mathrm{comp}}^{\mathrm{off}}V_{\mathrm{droop},0}
\end{aligned}
```

If optional reference inputs are not connected, use steady-state constants that
make the selected control errors zero. Omitted frequency inputs default to zero:

```math
\begin{aligned}
  f_0 &= 0 \\
  f_{\mathrm{ref},0} &= 0 \\
  e_{f,0} &= \text{deadband2}(f_{\mathrm{ref},0} - f_0, D_{f,\mathrm{bd1}}, D_{f,\mathrm{bd2}}) \\
  P_{\mathrm{freq},0} &= D_{\mathrm{dn}}\rho(e_{f,0}) - D_{\mathrm{up}}\rho(-e_{f,0}) \\
  V_{\mathrm{ref},0} &= V_{\mathrm{ctrl},0} \\
  Q_{\mathrm{ref},0} &= k_{\mathrm{base}}Q_{\mathrm{br},0} \\
  P_{\mathrm{plant},0}^{\mathrm{ref}} &= k_{\mathrm{base}}P_{\mathrm{br},0} - P_{\mathrm{freq},0}
\end{aligned}
```

Connected optional references use their supplied initial values and must satisfy the same steady-state residuals.

Initialize the measurement variables from the descriptor-form filter residuals:

```math
\begin{aligned}
  V_{\mathrm{meas},0} &= V_{\mathrm{ctrl},0} \\
  Q_{\mathrm{meas},0} &= k_{\mathrm{base}}Q_{\mathrm{br},0} \\
  P_{\mathrm{meas},0} &= k_{\mathrm{base}}P_{\mathrm{br},0}
\end{aligned}
```

Then evaluate the reactive-control algebraic chain:

```math
\begin{aligned}
  s_{\mathrm{frz},0} &= \text{above}(V_{\mathrm{reg},0}, V_{\mathrm{frz}}) \\
  e_{\mathrm{RQ},0}
    &= s_{\mathrm{ref}}\left(V_{\mathrm{ref},0} - V_{\mathrm{meas},0}\right)
       + s_{\mathrm{ref}}^{\mathrm{off}}\left(Q_{\mathrm{ref},0} - Q_{\mathrm{meas},0}\right) \\
  e_{\mathrm{RQ},0}^{\mathrm{db}}
    &= \text{deadband2}(e_{\mathrm{RQ},0}, D_{\mathrm{bd1}}, D_{\mathrm{bd2}}) \\
  e_{\mathrm{RQ},0}^{\mathrm{lim}}
    &= \text{clamp}(e_{\mathrm{RQ},0}^{\mathrm{db}}, e^{\min}, e^{\max})
\end{aligned}
```

Choose the initial reactive-command output and PI states as:

```math
\begin{aligned}
  Q_{\mathrm{ext},0} &= \text{clamp}(k_{\mathrm{base}}Q_{\mathrm{br},0}, Q^{\min}, Q^{\max}) \\
  Q_{\mathrm{PI},0} &= Q_{\mathrm{ext},0} \\
  x_{\mathrm{Qext},0} &= Q_{\mathrm{PI},0} \\
  x_{Q,0} &= Q_{\mathrm{ext},0} - K_{\mathrm{p}}e_{\mathrm{RQ},0}^{\mathrm{lim}}
\end{aligned}
```

This satisfies the lead-lag algebraic output and produces a zero lead-lag
state derivative. The reactive PI state derivative must also be zero after
antiwindup is applied.

If the `qext` output signal already carries a seeded value, REPCA uses that
value for $Q_{\mathrm{ext},0}$ and back-solves the same PI states from it.

Evaluate the active-power control chain:

```math
\begin{aligned}
  e_{P,0}
    &= P_{\mathrm{plant},0}^{\mathrm{ref}}
       - P_{\mathrm{meas},0}
       + P_{\mathrm{freq},0} \\
  e_{P,0}^{\mathrm{lim}} &= \text{clamp}(e_{P,0}, e_P^{\min}, e_P^{\max})
\end{aligned}
```

Choose:

```math
\begin{aligned}
  P_{\mathrm{ref},0} &= \text{clamp}(k_{\mathrm{base}}P_{\mathrm{br},0}, P^{\min}, P^{\max}) \\
  P_{\mathrm{PI},0} &= P_{\mathrm{ref},0} \\
  x_{P,0} &= P_{\mathrm{ref},0} - K_{\mathrm{pg}}e_{P,0}^{\mathrm{lim}} \\
  P_{\mathrm{ext},0} &= s_{\mathrm{freq}}P_{\mathrm{ref},0}
\end{aligned}
```

If the `pext` output signal already carries a seeded value and
$s_{\mathrm{freq}} \ne 0$, REPCA uses
$P_{\mathrm{ref},0}=P_{\mathrm{ext},0}/s_{\mathrm{freq}}$.

The initialized derivative vector is zero. Initialization rejects the case if
the reactive-power or active-power PI antiwindup rate is nonzero.

## Model Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`qext`          | [p.u.] | Reactive-power command output       | Component base; sent to downstream electrical control
`pext`          | [p.u.] | Active-power command output         | Component base; sent to downstream electrical control when `Freqflag = 1`
`vmeas`         | [p.u.] | Filtered regulated voltage          |
`qmeas`         | [p.u.] | Filtered reactive-power signal      |
`pmeas`         | [p.u.] | Filtered active-power signal        |
`pref`          | [p.u.] | Active-power command lag state      |
`vctrl`         | [p.u.] | Selected voltage-measurement input  |
`sfrz`          | [binary] | Reactive-PI voltage-enable indicator |
`qpi`           | [p.u.] | Reactive PI output                  |
`pfreq`         | [p.u.] | Computed frequency droop active-power correction |
`ppi`           | [p.u.] | Active-power PI output              |
