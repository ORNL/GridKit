# **Hydro Turbine-Governor Model (HYGOV)**

HYGOV is a hydro turbine-governor model with a temporary-droop governor, gate
servo, and nonlinear single-penstock water-column turbine. In GridKit it is
represented as a governor model that reads machine speed deviation and supplies
mechanical power to the machine.

Notes:
- HYGOV equations use the `Trate` base; the `pmech` signal uses `mva` and is converted at the signal boundary.
- GridKit requires `Trate > 0`; for PowerWorld `Trate = 0` records, set `Trate` to the connected machine `mva` before running.
- Aside from `db2`, model parameters are required input values; the table values are typical values, not parser defaults.
- HYGOV uses `db1` with CommonMath `deadband1`; HYGOVD-only `dbL`/`dbH` and nonzero `db2` are not modeled.
- PowerWorld fields `Ttur`, `Eps`, Kaplan blade-servo points, `Bmax`, and `Tblade` are not part of this implementation.
- Permanent-droop feedback uses desired gate $c$ (State 2), not gate $g$.

## Block Diagram

Standard HYGOV block diagram.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics/HYGOV_diagram.png">

  Figure 1: Governor HYGOV model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                          | Units    | JSON     | Description                                  | Typical Value | Note
--------------------------------|----------|----------|----------------------------------------------|---------------|------
$P^{\mathrm{rate}}$             | [MW]     | `Trate` | Component-base turbine rating                | 100.0         | Required positive value
$S^{\mathrm{pmech}}$            | [MVA]    | `mva`    | Mechanical-power signal base                 | 100.0         | Use connected machine `mva`
$R_{\mathrm{perm}}$             | [p.u.]   | `Rperm`  | Permanent droop                              | 0.04          | Source diagram label: `R`
$R_{\mathrm{temp}}$             | [p.u.]   | `Rtemp`  | Temporary droop                              | 0.3           | Source diagram label: `r`
$T_r$                           | [sec]    | `Tr`     | Temporary-droop reset time constant          | 5.0           |
$T_f$                           | [sec]    | `Tf`     | Governor error filter time constant          | 0.05          | State 1
$T_g$                           | [sec]    | `Tg`     | Gate servo time constant                     | 0.5           | State 3; if zero, $g$ is algebraic
$V_{\mathrm{elm}}$              | [p.u./s] | `Velm`   | Maximum desired-gate velocity magnitude      | 0.2           | Symmetric rate limit on State 2
$G^{\max}$                      | [p.u.]   | `Gmax`   | Maximum desired-gate position                | 1.0           |
$G^{\min}$                      | [p.u.]   | `Gmin`   | Minimum desired-gate position                | 0.0           |
$T_w$                           | [sec]    | `Tw`     | Water inertia time constant                  | 1.0           | State 4
$A_t$                           | [p.u.]   | `At`     | Turbine gain                                 | 1.2           |
$D_{\mathrm{turb}}$             | [p.u.]   | `Dturb`  | Turbine damping coefficient                  | 0.5           | Multiplied by speed deviation and gate
$q_{\mathrm{NL}}$               | [p.u.]   | `Qnl`    | No-load flow at nominal head                 | 0.05          |
$T_n$                           | [sec]    | `Tn`     | Speed lead-lag numerator time constant       | 0.0           | Source text labels this "lag"
$T_{np}$                        | [sec]    | `Tnp`    | Speed lead-lag denominator time constant     | 0.0           | Source text labels this "lead"
$D_{\omega}$                    | [p.u.]   | `db1`    | Type 1 speed deadband threshold              | 0.0           | Source data may provide this in Hz
$H_{\mathrm{dam}}$              | [p.u.]   | `Hdam`   | Head available at dam                        | 1.0           |

The nonlinear gate-to-power curve $N_{\mathrm{GV}}$ is represented by six source points:

Symbol                          | Units  | JSON     | Description                                  | Typical Value | Note
--------------------------------|--------|----------|----------------------------------------------|---------------|------
$G_V^{(k)}$                     | [p.u.] | `Gv0`-`Gv5`   | Gate point $k$ of the gain curve             | 0.0           | $k=0,\ldots,5$
$P_{\mathrm{GV}}^{(k)}$         | [p.u.] | `Pgv0`-`Pgv5` | Power point $k$ of the gain curve            | 0.0           | $k=0,\ldots,5$

If all six `Gv` and all six `Pgv` source values are zero, the source data is treated as
omitting the nonlinear curve. The effective curve used by the equations is the identity
curve with points `(0,0), (0.2,0.2), ..., (1,1)`.

### Parameter Validation

Invalid HYGOV parameter sets are rejected by the following checks. If source data preprocessing swaps gate limits, expands gate limits to include the initialized gate, adjusts small positive time constants, converts `db1` from Hz to p.u. speed, or replaces an all-zero source gate curve with the identity curve, apply these checks to the effective values used by the equations. Source data with nonzero `db2` is outside the supported equations below.

The required checks are:

```math
\begin{aligned}
  &P^{\mathrm{rate}} > 0 \\
  &S^{\mathrm{pmech}} > 0 \\
  &R_{\mathrm{perm}} > 0 \\
  &R_{\mathrm{temp}} > 0 \\
  &T_r > 0 \\
  &T_f > 0,\quad T_g \ge 0 \\
  &T_w > 0 \\
  &T_n \ge 0,\quad T_{np} \ge 0,\quad T_{np}=0 \Rightarrow T_n=0 \\
  &V_{\mathrm{elm}} > 0 \\
  &G^{\min} \le G^{\max} \\
  &A_t > 0 \\
  &D_{\mathrm{turb}} \ge 0 \\
  &q_{\mathrm{NL}} \ge 0 \\
  &D_{\omega} \ge 0 \\
  &H_{\mathrm{dam}} > 0 \\
  &G_V^{(0)} < G_V^{(1)} < \cdots < G_V^{(5)} \\
  &0 \le P_{\mathrm{GV}}^{(0)} \le P_{\mathrm{GV}}^{(1)} \le \cdots \le P_{\mathrm{GV}}^{(5)}
\end{aligned}
```

Initialization also requires $N_{\mathrm{GV}}^{-1}$ to be single-valued at the solved initial operating point.

### Model Derived Parameters

```math
\begin{aligned}
  k_n &=
    \begin{cases}
      T_n/T_{np} & T_{np} > 0 \\
      0 & T_{np} = 0
    \end{cases} \\
  S_{\mathrm{hygov}}^{\mathrm{base}}
    &= P^{\mathrm{rate}}
\end{aligned}
```

Here $S^{\mathrm{pmech}}$ is the `mva` parameter. The `pmech` signal uses
that base, while HYGOV turbine equations use the HYGOV component base:

```math
\begin{aligned}
  P^{\mathrm{hygov}}
    &= P^{\mathrm{pmech}}\dfrac{S^{\mathrm{pmech}}}{S_{\mathrm{hygov}}^{\mathrm{base}}} \\
  P^{\mathrm{pmech}}
    &= P^{\mathrm{hygov}}\dfrac{S_{\mathrm{hygov}}^{\mathrm{base}}}{S^{\mathrm{pmech}}}
\end{aligned}
```

The nonlinear gate-to-power function uses GridKit's smooth [Linear Segment](../../../../CommonMath.md#derived-functions) helper:

```math
\begin{aligned}
  N_{\mathrm{GV}}(x)
    &=
      P_{\mathrm{GV}}^{(0)}
      + \sum_{k\in\{0,\ldots,4\}}
        \text{linseg}\!\left(
          x;\,
          G_V^{(k)},\,
          G_V^{(k+1)},\,
          P_{\mathrm{GV}}^{(k+1)} - P_{\mathrm{GV}}^{(k)}
        \right)
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                         | Note
------------------------|--------|-------------------------------------|------
$x_n$                   | [p.u.] | Speed lead-lag denominator state    | Not circled in Fig. 1; realizes the `Tn`/`Tnp` block
$x_f$                   | [p.u.] | Governor error filter output        | State 1 in Fig. 1
$c$                     | [p.u.] | Desired-gate position               | State 2 in Fig. 1
$g$                     | [p.u.] | Gate position                       | State 3 in Fig. 1; algebraic when $T_g = 0$
$q$                     | [p.u.] | Turbine flow                        | State 4 in Fig. 1

#### Algebraic

Symbol                            | Units    | Description                         | Note
----------------------------------|----------|-------------------------------------|------
$\omega_{\mathrm{db}}$            | [p.u.]   | Type 1 deadbanded speed deviation   | Defined by CommonMath `deadband1`
$y_{\omega}$                      | [p.u.]   | Lead-lag-conditioned speed signal   | Output of the speed lead-lag
$e_f$                             | [p.u.]   | Governor error into the filter      | Reference path less conditioned speed and permanent-droop feedback
$f_c$                             | [p.u./s] | Desired-gate derivative target      | Before rate and position limits
$r_c$                             | [p.u./s] | Rate-limited desired-gate derivative target | Limited by $\pm V_{\mathrm{elm}}$
$P_{\mathrm{GV}}$                 | [p.u.]   | Nonlinear gate-to-power curve output | $N_{\mathrm{GV}}(g)$
$P_{\mathrm{GV}}^{\mathrm{safe}}$ | [p.u.]   | Safe gate-to-power value            | Lower bounded by 0.01
$H$                               | [p.u.]   | Turbine head                        | $(q/P_{\mathrm{GV}}^{\mathrm{safe}})^2$
$P_{\mathrm{mech}}$               | [p.u.]   | Mechanical power to generator       | `mva` signal; converted at the HYGOV boundary

### External Variables

#### Differential
None.

#### Algebraic

Symbol                          | Units  | Description                    | Note
--------------------------------|--------|--------------------------------|------
$\omega$                        | [p.u.] | Machine speed deviation        | Read from a machine model
$P_{\mathrm{ref}}$              | [p.u.] | Active-power/load reference    | HYGOV component base; external setpoint or constant parameter; source label: `Pref`
$P_{\mathrm{aux}}$              | [p.u.] | Auxiliary power input          | HYGOV component base; optional, defaults to zero; source label: `Paux`

## Model Equations

The desired-gate transfer block before limits is:

```math
\begin{aligned}
  \dfrac{c(s)}{x_f(s)}
    &= \dfrac{1+sT_r}{R_{\mathrm{temp}}T_rs}
\end{aligned}
```

The realization below forms the derivative target $f_c$, limits it to $r_c$, and applies anti-windup at the desired-gate limits.

### Differential Equations

```math
\begin{aligned}
  0 &= -T_{np}\dot x_n - x_n + \omega_{\mathrm{db}} \\
  0 &= -T_f\dot x_f - x_f + e_f \\
  0 &= -R_{\mathrm{temp}}T_r f_c + x_f + T_r\dot x_f \\
  0 &=
    -\dot c
    + \text{antiwindup}\!\left(
      c,
      r_c,
      G^{\min},
      G^{\max}
    \right) \\
  0 &= -T_g\dot g - g + c \\
  0 &= -T_w\dot q + H_{\mathrm{dam}} - H
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#anti-windup-indicator) target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -\omega_{\mathrm{db}}
       + \text{deadband1}\!\left(\omega, -D_{\omega}, D_{\omega}\right) \\
  0 &= -y_{\omega}
       + x_n
       + k_n\left(\omega_{\mathrm{db}} - x_n\right) \\
  0 &= -e_f + P_{\mathrm{ref}} + P_{\mathrm{aux}} - y_{\omega} - R_{\mathrm{perm}}c \\
  0 &= -r_c + \text{clamp}\!\left(f_c, -V_{\mathrm{elm}}, V_{\mathrm{elm}}\right) \\[0.5ex]
  0 &= -P_{\mathrm{GV}} + N_{\mathrm{GV}}(g) \\
  0 &= -P_{\mathrm{GV}}^{\mathrm{safe}} + \text{max}\!\left(P_{\mathrm{GV}}, 0.01\right) \\
  0 &= -H + \left(\dfrac{q}{P_{\mathrm{GV}}^{\mathrm{safe}}}\right)^2 \\
  0 &= -P_{\mathrm{mech}}^{\mathrm{pmech}}
       + \left(A_tH(q - q_{\mathrm{NL}}) - D_{\mathrm{turb}}\omega g\right)
         \dfrac{S_{\mathrm{hygov}}^{\mathrm{base}}}{S^{\mathrm{pmech}}}
\end{aligned}
```

The $P_{\mathrm{GV}}^{\mathrm{safe}}$ lower bound protects the head divider near a closed gate. A physically consistent operating point satisfies $P_{\mathrm{GV}}>0.01$ so the lower bound is inactive.

CommonMath defines the helper targets and smooth approximations for [max, clamp, deadband1, and linseg](../../../../CommonMath.md#derived-functions).

## Initialization

Initialization is performed by evaluating the steady-state residuals in dependency order. Let subscript $0$ denote initial values and set all internal derivatives to zero. If optional signals are not connected, use:

```math
\begin{aligned}
  \omega_0 &= 0 \\
  P_{\mathrm{aux},0} &= 0
\end{aligned}
```

The connected power-source model initializes mechanical power first; HYGOV reads
the `pmech` signal $P_{\mathrm{mech},0}^{\mathrm{pmech}}$ and converts it to
the component-base value $P_{\mathrm{mech},0}^{\mathrm{hygov}}$. At synchronous
speed the damping term vanishes:

```math
\begin{aligned}
  P_{\mathrm{mech},0}^{\mathrm{hygov}}
    &= P_{\mathrm{mech},0}^{\mathrm{pmech}}
       \dfrac{S^{\mathrm{pmech}}}{S_{\mathrm{hygov}}^{\mathrm{base}}} \\
  H_0 &= H_{\mathrm{dam}} \\
  q_0 &= q_{\mathrm{NL}} + \dfrac{P_{\mathrm{mech},0}^{\mathrm{hygov}}}{A_tH_0} \\
  P_{\mathrm{GV},0} &= \dfrac{q_0}{\sqrt{H_0}} \\
  P_{\mathrm{GV},0}^{\mathrm{safe}} &= \text{max}\!\left(P_{\mathrm{GV},0}, 0.01\right) \\
  g_0 &= N_{\mathrm{GV}}^{-1}\!\left(P_{\mathrm{GV},0}\right) \\
  c_0 &= g_0
\end{aligned}
```

Then evaluate the governor chain. With $\dot c_0=0$ and $\dot x_{f,0}=0$, the temporary-droop block requires $x_{f,0}=0$:

```math
\begin{aligned}
  \omega_{\mathrm{db},0} &= \text{deadband1}\!\left(\omega_0, -D_{\omega}, D_{\omega}\right) \\
  x_{n,0} &= \omega_{\mathrm{db},0} \\
  y_{\omega,0} &= \omega_{\mathrm{db},0} \\
  x_{f,0} &= 0 \\
  e_{f,0} &= x_{f,0} \\
  f_{c,0} &= 0 \\
  r_{c,0} &= 0 \\
  P_{\mathrm{ref},0} &= e_{f,0} - P_{\mathrm{aux},0} + y_{\omega,0} + R_{\mathrm{perm}}c_0
\end{aligned}
```

For an unsaturated zero-derivative start, require $G^{\min} \le c_0 \le G^{\max}$, $P_{\mathrm{GV},0}>0.01$, and a single-valued inverse $N_{\mathrm{GV}}^{-1}(P_{\mathrm{GV},0})$. Starts that bind a gate limit or cannot be served at the available head are outside these closed-form initialization formulas.

## Model Outputs

Output         | Units  | Description                         | Note
---------------|--------|-------------------------------------|------
`pmech`        | [p.u.] | Mechanical power output             | `mva` signal; oriented as mechanical input to the machine
`leadlag`      | [p.u.] | Speed lead-lag denominator state    | Not circled in Fig. 1
`filter`       | [p.u.] | Governor error filter output        | State 1
`desiredgate`  | [p.u.] | Desired-gate position               | State 2
`gate`         | [p.u.] | Gate position                       | State 3
`flow`         | [p.u.] | Turbine flow                        | State 4
`head`         | [p.u.] | Turbine head                        | $(q/P_{\mathrm{GV}}^{\mathrm{safe}})^2$
`pgv`          | [p.u.] | Nonlinear gate-to-power curve output |
