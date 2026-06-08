# **Renewable Energy Generator/Converter Model (REGCA)**

REGCA is a first-generation WECC renewable generator/converter model for
inverter-coupled resources. In GridKit it is represented as a controlled
current source at the network interface.

Notes:
- Internal currents are on converter base.
- Bus injections and exported branch powers are on system base.
- LVACM uses $V_T$; LVPL uses $V_M$.
- HVRCM is represented by algebraic current $I_{\mathrm{q}}^{\mathrm{extra}}$.

## Block Diagram

Standard REGCA converter-interface model.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_REGCA_Diagram.png">

  Figure 1: Generator/Converter REGCA model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                           | Units    | JSON     | Description                                           | Typical Value | Note
---------------------------------|----------|----------|-------------------------------------------------------|---------------|------
$P_0$                            | [p.u.]   | `P0`     | Initial active power injection on system base         | 1.0           | Required initialization source
$Q_0$                            | [p.u.]   | `Q0`     | Initial reactive power injection on system base       | 0.0           | Required initialization source
$S^{\mathrm{base}}$              | [MVA]    | `mva`    | REGCA model power base                                | 100.0         |
$T_{\mathrm{g}}$                 | [sec]    | `Tg`     | Converter current-control lag time constant           | 0.02          | Block name: `Tg`
$T_M$                            | [sec]    | `TM`     | Terminal voltage sensor time constant                 | 0.02          | Block name: `Tfltr`
$R_{\mathrm{q}}^{\max}$           | [p.u./s] | `Rqmax`  | Reactive-current recovery positive rate limit         | 999.0         | Block name: `Iqrmax`; GridKit requires $R_{\mathrm{q}}^{\max} > 0$
$R_{\mathrm{q}}^{\min}$           | [p.u./s] | `Rqmin`  | Reactive-current recovery negative rate limit         | -999.0        | Block name: `Iqrmin`; GridKit requires $R_{\mathrm{q}}^{\min} < 0$
$R_{\mathrm{p}}^{\max}$           | [p.u./s] | `Rpmax`  | Active-current magnitude recovery rate limit          | 999.0         | Block name: `rrpwr`
$s_L$                            | [binary] | `sL`     | LVPL switch                                           | 1             | Block name: `LPVLSW`
$I_{L1}$                         | [p.u.]   | `IL1`    | LVPL upper-current ceiling                            | 1.1           | Block name: `LVPL1`
$V_{L0}$                         | [p.u.]   | `VL0`    | LVPL zero-crossing voltage                            | 0.4           | Block name: `zerox`
$V_{L1}$                         | [p.u.]   | `VL1`    | LVPL upper breakpoint voltage                         | 0.9           | Block name: `brkpt`
$V_{A0}$                         | [p.u.]   | `VA0`    | LVACM lower breakpoint voltage                        | 0.4           | Block name: `LVPnt0`
$V_{A1}$                         | [p.u.]   | `VA1`    | LVACM upper breakpoint voltage                        | 0.9           | Block name: `LVPnt1`
$V_{\mathrm{hv}}^{\max}$          | [p.u.]   | `Vhvmax` | Terminal-voltage ceiling for HV reactive management   | 1.2           | Block name: `VLim`

PowerWorld fields `Qmin`, `Khv`, and `Xe` are not parameters of this GridKit
REGCA implementation.

### Parameter Validation

Implementations should reject or report invalid parameter sets:

```math
\begin{aligned}
  S^{\mathrm{base}} &> 0 &
  T_{\mathrm{g}} &> 0 &
  T_M &> 0 \\
  R_{\mathrm{p}}^{\max} &> 0 &
  R_{\mathrm{q}}^{\min} &< 0 < R_{\mathrm{q}}^{\max} &
  s_L &\in \{0,1\} \\
  I_{L1} &\ge 0 &
  0 &\le V_{L0} < V_{L1} &
  0 &\le V_{A0} < V_{A1} \\
  V_{\mathrm{hv}}^{\max} &> 0
\end{aligned}
```

### Model Derived Parameters

The active-current bounds use finite surrogate $M_{\mathrm{p}}$ for inactive
$\pm\infty$ limits:

```math
M_{\mathrm{p}} = 100 R_{\mathrm{p}}^{\max}
```

$M_{\mathrm{p}}$ is not a physical REGCA parameter.

## Model Variables

### Internal Variables

#### Differential

Symbol                | Units  | Description               | Note
----------------------|--------|---------------------------|------
$V_M$                 | [p.u.] | Filtered terminal voltage | State 3 in Fig. 1
$I_{\mathrm{q}}$      | [p.u.] | Reactive-current state    | State 1 in Fig. 1 before the `-1` block; converter base
$I_{\mathrm{p}}$      | [p.u.] | Active-current state      | State 2 in Fig. 1; converter base

#### Algebraic

Symbol                     | Units  | Description                                                 | Note
---------------------------|--------|-------------------------------------------------------------|------
$V_T$                      | [p.u.] | Terminal voltage magnitude                                  |
$I_{\mathrm{i}}$           | [p.u.] | Injected current, imaginary component on network reference frame | Converter base
$I_{\mathrm{q}}^{\mathrm{extra}}$ | [p.u.] | Extra inductive current from high-voltage reactive current management | Converter base
$I_L$                      | [p.u.] | LVPL upper-limit current curve                              | Function of $V_M$
$I_{\mathrm{r}}$           | [p.u.] | Injected current, real component on network reference frame | Converter base
$\ell_{\mathrm{p}}$        | [p.u./s] | Smooth active-current lower rate bound                    | Equivalent to diagram `Rdown`
$u_{\mathrm{p}}$           | [p.u./s] | Smooth active-current upper rate bound                    | Effective `Rup`; includes LVPL anti-windup when $s_L=1$
$P_{\mathrm{br}}$          | [p.u.] | Active-power output                                         | System base
$Q_{\mathrm{br}}$          | [p.u.] | Reactive-power output                                       | System base

### External Variables

#### Differential
None.

#### Algebraic

Symbol                          | Units  | Description                                                      | Note
--------------------------------|--------|------------------------------------------------------------------|------
$V_{\mathrm{r}}$                | [p.u.] | Terminal voltage, real component on network reference frame      | Owned by bus object
$V_{\mathrm{i}}$                | [p.u.] | Terminal voltage, imaginary component on network reference frame | Owned by bus object
$I_{\mathrm{p}}^{\mathrm{cmd}}$ | [p.u.] | Active-current command in the terminal-voltage reference frame   | Converter base; owned by REEC, constant if no REEC is connected
$I_{\mathrm{q}}^{\mathrm{cmd}}$ | [p.u.] | Reactive-current command in the terminal-voltage reference frame | Converter base; owned by REEC, constant if no REEC is connected

REGCA initializes its current commands from the required `P0`/`Q0` power-flow
injection and writes those resolved commands to any attached `ipcmd`/`iqcmd`
ports. If no controller is connected, the resolved initialization commands are
held constant during residual evaluation.

## Model Equations

### Differential Equations

The state equations use CommonMath helper notation. The $I_{\mathrm{q}}$
limiter branch is selected by the initialized reactive-current command.

```math
\begin{aligned}
  0 &= -T_M \dot V_M - V_M + V_T \\
  0 &= -T_{\mathrm{g}}\dot I_{\mathrm{q}} +
       (I_{\mathrm{q}}^{\mathrm{cmd}} - I_{\mathrm{q}}) +
    \begin{cases}
      -\rho(I_{\mathrm{q}}^{\mathrm{cmd}} - I_{\mathrm{q}} - T_{\mathrm{g}}R_{\mathrm{q}}^{\max})
        & I_{\mathrm{q},0}^{\mathrm{cmd}} > 0 \\
      \rho(T_{\mathrm{g}}R_{\mathrm{q}}^{\min} - (I_{\mathrm{q}}^{\mathrm{cmd}} - I_{\mathrm{q}}))
        & I_{\mathrm{q},0}^{\mathrm{cmd}} \le 0
    \end{cases} \\
  0 &= -T_{\mathrm{g}}\dot I_{\mathrm{p}}
       + T_{\mathrm{g}}\ell_{\mathrm{p}}
       + \rho(I_{\mathrm{p}}^{\mathrm{cmd}} - I_{\mathrm{p}} - T_{\mathrm{g}}\ell_{\mathrm{p}})
       - \rho(I_{\mathrm{p}}^{\mathrm{cmd}} - I_{\mathrm{p}} - T_{\mathrm{g}}u_{\mathrm{p}})
\end{aligned}
```

CommonMath defines the [linear segment, ramp, and sigmoid helpers](../../../../CommonMath.md#derived-functions).

### Algebraic Equations

The algebraic equations are:

```math
\begin{aligned}
  0 &= -V_T^2 + V_{\mathrm{r}}^2 + V_{\mathrm{i}}^2 \\
  0 &= -V_T I_{\mathrm{i}}
       - (I_{\mathrm{q}} - I_{\mathrm{q}}^{\mathrm{extra}})V_{\mathrm{r}}
       + V_{\mathrm{i}}I_{\mathrm{p}}\text{linseg}(V_T;\ V_{A0},\ V_{A1},\ 1) \\
  0 &= -I_{\mathrm{q}}^{\mathrm{extra}}
       + \rho(V_T - V_{\mathrm{hv}}^{\max}) \\
  0 &= -I_L
       + \text{linseg}(V_M;\ V_{L0},\ V_{L1},\ I_{L1}) \\
  0 &= -V_T I_{\mathrm{r}}
       + (I_{\mathrm{q}} - I_{\mathrm{q}}^{\mathrm{extra}})V_{\mathrm{i}}
       + V_{\mathrm{r}}I_{\mathrm{p}}\text{linseg}(V_T;\ V_{A0},\ V_{A1},\ 1) \\
  0 &= -\ell_{\mathrm{p}}
       - R_{\mathrm{p}}^{\max}
       - (M_{\mathrm{p}} - R_{\mathrm{p}}^{\max})\sigma(I_{\mathrm{p}}) \\
  0 &= -u_{\mathrm{p}}
       + M_{\mathrm{p}}(1-\sigma(I_{\mathrm{p}}))
       + R_{\mathrm{p}}^{\max}\sigma(I_{\mathrm{p}})
         \left(1 - s_L + s_L\sigma(I_L - I_{\mathrm{p}})\right) \\
  0 &= -P_{\mathrm{br}}
       + \dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}}
         \left(V_{\mathrm{r}} I_{\mathrm{r}} + V_{\mathrm{i}} I_{\mathrm{i}}\right) \\
  0 &= -Q_{\mathrm{br}}
       + \dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}}
         \left(V_{\mathrm{i}} I_{\mathrm{r}} - V_{\mathrm{r}} I_{\mathrm{i}}\right)
\end{aligned}
```

The $V_T$ residual is kept in squared form for smoothness at the origin.

## Network Interface

The bus receives system-base current injections converted from converter-base
REGCA currents:

```math
\begin{aligned}
  I_{\mathrm{r}}^{\mathrm{inj}} &:= I_{\mathrm{r}}\dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}} \\
  I_{\mathrm{i}}^{\mathrm{inj}} &:= I_{\mathrm{i}}\dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}}
\end{aligned}
```

Positive current injection is into the bus.

## Initialization

Given initialized bus voltage $V_{\mathrm{r}}, V_{\mathrm{i}}$, compute the
steady-state initial values:

```math
\begin{aligned}
  V_T                         &= \sqrt{V_{\mathrm{r}}^2 + V_{\mathrm{i}}^2} \\
  V_{M,0}                      &= V_T \\
  I_{L,0}                      &= \text{linseg}(V_T;\ V_{L0},\ V_{L1},\ I_{L1}) \\
  I_{\mathrm{p},0}             &= V_T^{-1}P_0
                                  \dfrac{S^{\mathrm{sys}}}{S^{\mathrm{base}}}
                                  \left[\text{linseg}(V_T;\ V_{A0},\ V_{A1},\ 1)\right]^{-1} \\
  I_{\mathrm{q},0}             &= V_T^{-1}Q_0\dfrac{S^{\mathrm{sys}}}{S^{\mathrm{base}}} \\
  \ell_{\mathrm{p},0}          &= -R_{\mathrm{p}}^{\max}
       - (M_{\mathrm{p}} - R_{\mathrm{p}}^{\max})\sigma(I_{\mathrm{p},0}) \\
  u_{\mathrm{p},0}             &= M_{\mathrm{p}}(1-\sigma(I_{\mathrm{p},0}))
       + R_{\mathrm{p}}^{\max}\sigma(I_{\mathrm{p},0})
         \left(1 - s_L + s_L\sigma(I_{L,0} - I_{\mathrm{p},0})\right) \\
  I_{\mathrm{q},0}^{\mathrm{extra}} &= 0 \\
  I_{\mathrm{i},0}             &= V_T^{-1}\left[
       -I_{\mathrm{q},0}V_{\mathrm{r}}
       + V_{\mathrm{i}}I_{\mathrm{p},0}\text{linseg}(V_T;\ V_{A0},\ V_{A1},\ 1)
     \right] \\
  I_{\mathrm{r},0}             &= V_T^{-1}\left[
       I_{\mathrm{q},0}V_{\mathrm{i}}
       + V_{\mathrm{r}}I_{\mathrm{p},0}\text{linseg}(V_T;\ V_{A0},\ V_{A1},\ 1)
     \right] \\
  P_{\mathrm{br},0}            &= \dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}}
                                  \left(V_{\mathrm{r}} I_{\mathrm{r},0} + V_{\mathrm{i}} I_{\mathrm{i},0}\right) \\
  Q_{\mathrm{br},0}            &= \dfrac{S^{\mathrm{base}}}{S^{\mathrm{sys}}}
                                  \left(V_{\mathrm{i}} I_{\mathrm{r},0} - V_{\mathrm{r}} I_{\mathrm{i},0}\right)
\end{aligned}
```

REGCA writes the resolved initial commands to attached `ipcmd` and `iqcmd`
ports. Initialization rejects nonpositive $V_T$, terminal voltage at or above
$V_{\mathrm{hv}}^{\max}$, and nonzero $P_0$ with zero LVACM gain. All internal
derivatives initialize to zero.

## Model Outputs

Output          | Units    | Description                          | Note
----------------|----------|--------------------------------------|------
`ir`            | [p.u.]   | Real current injection               | Converter base; exported through `ibranchr` when assigned
`ii`            | [p.u.]   | Imaginary current injection          | Converter base; exported through `ibranchi` when assigned
`p`             | [p.u.]   | Active-power output                  | System base; exported through `pbranch` when assigned
`q`             | [p.u.]   | Reactive-power output                | System base; exported through `qbranch` when assigned
`vt`            | [p.u.]   | Terminal voltage magnitude           |
`vm`            | [p.u.]   | Filtered terminal voltage            |
`ip`            | [p.u.]   | Active-current state                 | Converter base
`iq`            | [p.u.]   | Reactive-current state               | Converter base
`iqextra`       | [p.u.]   | HVRCM extra reactive current         | Converter base
`il`            | [p.u.]   | LVPL upper-limit current curve       |
`lp`            | [p.u./s] | Active-current lower rate bound      |
`up`            | [p.u./s] | Active-current upper rate bound      |
