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

The IEEET1 differential equations, as derived from the model diagram.  By defining an auxillary valve function $f=-V_{R}+K_{a}V_{tr}$ one can write the piecewise definition of $\dot V_R$ compactly.
```math
\begin{aligned}
   \dot{V}_{ts}   &= \dfrac{1}{T_R}(E_C-V_{ts}) \\
   \dot{V}_{R}    &= \dfrac{1}{T_A}
   \begin{cases}
      -V_{R} + K_{a}V_{tr}
         &  \text{if } (V_{rmin} < V_R < V_{rmax}) & \lor \\
         &  \quad (V_R \leq V_{rmin} \land f>0)    & \lor \\
         &  \quad(V_R \geq V_{rmax} \land f<0)            \\
      0  &  \text{else } \\
   \end{cases} \\
   \dot{E}_{fd}'  &= \dfrac{1}{T_E}(V_R-V_E-K_E E_{fd}') \\
   \dot{V}_{fx}   &= \dfrac{1}{T_F}V_f \\
\end{aligned}
```

#### Smooth Piecewise Approximation (Differential) 

The domain of the state variable $V_{R}\in(V_{rmin}, V_{rmax})$ is enforced
through the piece-wise definition above. This may need to be expressed as a
smooth approximation (smooth indicator $\phi$) expressed generically as follows.
```math
\begin{aligned}
   f(V_R)     
            &= 
            \dfrac{1}{T_A}
            \left[
               -V_{R}+K_{a}V_{tr}
            \right] \\
   \dot{V}_R &= 
            \phi(V_R, f) \cdot f(V_R)
\end{aligned}
```

The indicator function $\phi$ can be defined in terms of a scaled activation function ($\sigma$, sigmoid) and the $V_R$ limits $(V_{rmin}, V_{rmax})$.
```math
\begin{aligned}
   \phi_L(V_R,f)&= \sigma(V_{rmin}-V_R)\sigma(-f)            \\
   \phi_U(V_R,f)&= \sigma(V_R-V_{rmax})\sigma(f)             \\
   \phi(V_R,f)  &= \left[1-\phi_L\right]\left[1-\phi_U\right]\\
\end{aligned}
```

The scale of the sigmoid function ($\alpha$ on the order of $10^3$) should be chosen so that for all practical parameters of the IEEET1 model, the sigmoid acts as a step function.
```math
\begin{aligned}
   \sigma(x) = 
      \dfrac{1}{1+\exp(-\alpha x)}
\end{aligned}
```

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

At steady state we assume that $V_R$ is at or within its limits, 
and that the exciter is not saturated. The field $E_{fd}$ can be 
obtained through the steady-state conditions of the machine. 
We also assume for the moment that the stabilizer and over/under 
excitation limiters are non-existant. As of now, we assume there is no compensation impedence and that $E_C$ is simply the terminal voltage magnitude.
```math
\begin{aligned}
    E_{fd}' &= E_{fd}  \\
    V_R     &= K_E E_{fd} \\
    V_{fx}  &= \dfrac{K_F}{T_F} E_{fd} \\
    V_{tr}  &= \dfrac{K_E}{K_{a}} E_{fd} \\
    E_C     &= V_{term}\\
    V_{ts}  &= V_{term}\\
    V_{ref} &= E_c + V_{tr} \\
    V_{f}   &= 0 \\
    V_{E}   &= 0 \\
    k_{sat} &= 0 \\
\end{aligned}
```
