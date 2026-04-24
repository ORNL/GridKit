# **IEEE Type 1 Excitation System Model (IEEET1)**

## Block Diagram

Standard model of the IEEET1 Exciter.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_IEEET1_Diagram.png">
   
   
  Figure 1: Exciter IEEET1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters


Symbol      | Units  | Description                          | Typical Value | Note
------------|--------|--------------------------------------|---------| ------
$T_R$       | [sec]  | Time constant for voltage sensing    | 0       |
$K_a$       | [p.u.] | Coefficient for voltage regulation   | 50      |
$T_a$       | [sec]  | Time constant for voltage regulation | 0.04    |
$K_e$       | [p.u.] | Coefficient for excitation system    | -0.06   |
$T_e$       | [sec]  | Time constant for excitation system  | 0.6     | 
$K_f$       | [p.u.] | Coefficient for feedback             | 0.09    | 
$T_f$       | [sec]  | Time constant for feedback           | 1.46    | 
$V_{rmin}$  | [p.u.] | Lower limit to voltage regulation    | -1      | 
$V_{rmax}$  | [p.u.] | Upper limit to voltage regulation    | 1       | 
$E_1$       | [p.u.] | Saturation Parameter                 | 2.8     | 
$E_2$       | [p.u.] | Saturation Parameter                 | 3.73    | 
$S_{e1}$    | [p.u.] | Saturation Parameter                 | 0.04    | 
$S_{e2}$    | [p.u.] | Saturation Parameter                 | 0.33    | 
$I_{spdlm}$ | [binary] | Speed Limit flag indicator         | 0       | 

### Model Derived Parameters

The relationship of the derived parameters is defined by the following quadratic model. The parameters are chosen so that the quadratic model represents
the expected saturation near the operating region.
``` math
\begin{aligned}
  S_{e1} &= S_B(E_1-S_A)^2 \\
  S_{e2} &= S_B(E_2-S_A)^2 \\
\end{aligned}
```

Generally, this system has two solutions. The non-extraneous solution is as follows.
``` math
\begin{aligned}
  C &=  \sqrt{
   \dfrac
   {S_{e2}}
   {S_{e1}}
  } 
  \\
  S_A &= 
   \dfrac
   {C E_1 - E_2}
   {C - 1} 
  \\
  S_B &= 
   \dfrac
   {S_{e1}}
   {(E_1-S_A)^2}
\end{aligned}
```

## Model Variables 

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-------
$V_{ts}$  | [p.u.] | Sensed terminal voltage           |
$V_{R}$   | [p.u.] | Voltage regulator                 | 
$E_{fd}'$ | [p.u.] | Field-current pre-speed multiplier| 
$V_{fx}$  | [p.u.] | Exciter feedback internal state   | 


#### Algebraic


Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$V_{tr}$        | [p.u.] | Terminal Voltage Error            |
$V_{f}$         | [p.u.] | Feedback Voltage                  |
$V_{E}$         | [p.u.] | Excitation control voltage        |
$E_{fd}$        | [p.u.] | Field winding voltage             |
$k_{sat}$       | [p.u.] | Saturation variable               |

### External Variables

#### Differential
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$\omega$  | [p.u.] | Machine Speed Deviation                   | Read from a Machine Model

#### Algebraic


Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$E_{C}$         | [p.u.] | Compensated machine terminal voltage magnitude | 
$V_{ref}$       | [p.u.] | Reference terminal voltage                     |
$V_{UEL}$       | [p.u.] | Input from under excitation limiter            |
$V_{OEL}$       | [p.u.] | Input from over excitation limiter             |
$V_{S}$         | [p.u.] | Input from stabilizer controller               |


## Model Equations

### Differential Equations

The IEEET1 differential equations, as derived from the model diagram. Define the pre-limit derivative of $V_R$

```math
f = \dfrac{1}{T_A}\left[-V_R + K_a V_{tr}\right]
```

so that $\dot V_R$ can be written in piecewise form compactly.
```math
\begin{aligned}
   \dot{V}_{ts}   &= \dfrac{1}{T_R}(E_C-V_{ts}) \\
   \dot{V}_{R}    &=
   \begin{cases}
      f
         &  \text{if } (V_{rmin} < V_R < V_{rmax}) & \lor \\
         &  \quad (V_R \leq V_{rmin} \land f>0)    & \lor \\
         &  \quad(V_R \geq V_{rmax} \land f<0)            \\
      0  &  \text{else } \\
   \end{cases} \\
   \dot{E}_{fd}'  &= \dfrac{1}{T_E}(V_R-V_E-K_E E_{fd}') \\
   \dot{V}_{fx}   &= \dfrac{1}{T_F}V_f \\
\end{aligned}
```

In simulation the piecewise form above is replaced with a smooth approximation where $\phi$ is GridKit's smooth anti-windup indicator. See [CommonMath: Anti-Windup Indicator](../../../../CommonMath.md#anti-windup-indicator) for its definition, behavior, and design rationale.

### Algebraic Equations

The algebraic equations of the exciter.
```math
\begin{aligned}
   V_{tr} &= V_{ref} +V_{UEL} + V_{OEL} + V_S - V_f- V_{ts}\\
    V_{f} &= \dfrac{E_{fd}' K_F}{T_F} - V_{fx}\\
    V_{E} &= k_{sat}\cdot E_{fd}' \\
    E_{fd}&= \begin{cases}
        (1+\omega)E_{fd}'           &  \text{if } I_{spdlm}=1\\
         E_{fd}' &  \text{else} \\
   \end{cases}\\
    k_{sat}&= \begin{cases}
        S_B(E_{fd}' -S_A)^2        &  \text{if } E_{fd}' >S_A\\
        0  &  \text{else } \\
   \end{cases} \\
\end{aligned}
```
#### Smooth Piecewise Approximation (Algebraic) 

For the algebraic piecewise functions (non-flags), this implementation is straightforward when the approximation above is used.
```math
\begin{aligned}
    E_{fd}
    &=(1 + \omega I_{spdlm})E_{fd}' \\ 
    k_{sat}
    &=S_B\left[(E_{fd}' -S_A) \cdot \sigma (E_{fd}' -S_A)\right]^2   
\end{aligned}
```

The approximation approaches an exact solution as $\alpha\to\infty$.

## Initialization

The machine initializes $E_{fd}$ first. IEEET1 reads that value as $E_{fd,0}$, along with any attached $\omega$ and $V_S$, and solves the steady-state algebraic chain so all residuals vanish with $\dot y = 0$. There is no compensation impedance, so $E_C$ is taken as the terminal voltage magnitude. Saturation and the speed-limit flag are included directly; $V_{ref}$ is set to close the $V_{tr}$ equation with the current auxiliary inputs.

```math
\begin{aligned}
   E_C      &= \sqrt{V_r^2 + V_i^2} \\
   E_{fd}'  &= \dfrac{E_{fd,0}}{1 + I_{spdlm}\,\omega} \\
   k_{sat}  &= S_B\big[(E_{fd}' - S_A)\,\sigma(E_{fd}' - S_A)\big]^2 \\
   V_E      &= k_{sat}\, E_{fd}' \\
   V_R      &= K_E\, E_{fd}' + V_E \\
   V_{tr}   &= \dfrac{V_R}{K_a} \\
   V_{fx}   &= \dfrac{K_F}{T_F}\, E_{fd}' \\
   V_{ts}   &= E_C \\
   V_f      &= 0 \\
   V_{ref}  &= E_C + V_{tr} - V_{UEL} - V_{OEL} - V_S
\end{aligned}
```

All internal derivatives initialize to zero.
## Model Outputs

The field voltage, $E_{fd}$, is an internal model variable.

The magnetic saturation coefficient $k_{sat}$ is calculated from $E_{fd}$ using the the smooth piecewise version (above).
