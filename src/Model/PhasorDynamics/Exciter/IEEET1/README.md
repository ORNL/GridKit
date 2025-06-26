# **IEEE Type 1 Excitation System Model (IEEET1)**

## Block Diagram

Standard model of the IEEET1 Exciter.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_IEEET1_Diagram.png">
   
   
  Figure 1: Exciter IEEET1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters


Symbol      | Units  | Description                       | Typical Value | Note
------------|--------|-----------------------------------|---------------| ------
$T_R$       | [sec]  | Time constant for voltage sensing                 | 0       |
$K_a$       | [p.u.] | Coefficient for voltage regulation                | 50      |
$T_a$       | [sec]  | Time constant for voltage regulation              | 0.04    |
$K_e$       | [p.u.] | Coefficient for excitation system                 | -0.06   |
$T_e$       | [sec]  | Time constant for excitation system               | 0.6     | 
$K_f$       | [p.u.] | Tand time constant                                | 0.09    | 
$T_f$       | [sec]  | Coefficient for feedback                          | 1.46    | 
$V_{rmin}$  | [p.u.] | Lower limit to voltage regulation                 | -1      | 
$V_{rmax}$  | [p.u.] | Upper limit to voltage regulation                 | 1       | 
$E_1$       | [p.u.] | Saturation Parameter                              | 2.8     | 
$E_2$       | [p.u.] | Saturation Parameter                              | 3.73    | 
$S_{e1}$    | [p.u.] | Saturation Parameter                              | 0.04    | 
$S_{e2}$    | [p.u.] | Saturation Parameter                              | 0.33    | 
$I_{spdlm}$ | [binary] | Speed Limit flag indicator                        | 0       | 


## Model Variables 

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-------
$V_{ts}$  | [p.u.] | Sensed terminal voltage           |
$V_{R}$   | [p.u.] | Voltage regulator                 | 
$E_{fd}'$ | [p.u.] | Field-current pre-speed multiplier| 
$V_{fx}$ | [p.u.] | Exciter feedback internal state   | 


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
$\Delta\omega$  | [p.u.] | Speed Deviation                   | Read from a Machine Model

#### Algebraic


Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$E_{C}$         | [p.u.] | Compensated machine terminal voltage magnitude  | 
$V_{ref}$       | [p.u.] | Reference terminal voltage                   |
$V_{UEL}$       | [p.u.] | Input from under excitation limiter                   |
$V_{OEL}$       | [p.u.] | Input from over excitation limiter                   |
$V_{S}$         | [p.u.] | Input from stabilizer controller                  |


## Model Equations

### Differential Equations

The IEEET1 differential equations, as derived from the model diagram.  By defining an auxillary valve function $f=-V_{R}+K_{a}V_{tr}$ one can write the piecewise definition of $\dot V_R$ compactly.
```math
\begin{aligned}
   \dot{V}_{ts}   &= \dfrac{1}{T_R}(E_C-V_{ts}) \\
   \dot{V}_{R}    &= 
      \dfrac{1}{T_A}
   \begin{cases}
      -V_{R}+K_{a}V_{tr}
         &  \text{if } (V_{rmin} < V_R < V_{rmax}) & \lor \\
         &  \quad (V_R \leq V_{rmin} \land f>0)    & \lor \\
         &  \quad(V_R \geq V_{rmax} \land f<0)            \\
      0  
         &  \text{else } \\
   \end{cases} \\
   
   \dot{E}_{fd}'  &= \dfrac{1}{T_E}(V_R-V_E-K_E E_{fd}') \\
   \dot{V}_{fx}   &= \dfrac{1}{T_F}V_f \\
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
        (1+\Delta \omega)E_{fd}'           &  \text{if } I_{spdlm}=1\\
         E_{fd}' &  \text{else} \\
   \end{cases}\\
    k_{sat}&= \begin{cases}
        S_B(E_{fd}' -S_A)^2        &  \text{if } E_{fd}' >S_A\\
        0  &  \text{else } \\
   \end{cases} \\
\end{aligned}
```


#### Smooth Piecewise Approximation 

The domain of the state variable $V_{R}\in(V_{rmin}, V_{rmax})$ is enforced
through the piece-wise definition above. This may need to be expressed as a
smooth approximation (smooth indicator $\phi$) expressed generically as follows.
```math
\begin{aligned}
   f(V_R)     
            &= 
            \phi(V_R)\cdot \dfrac{1}{T_A}
            \left[
               -V_{R}+K_{a}V_{tr}
            \right] \\
   \dot{V}_R &= 
            \phi(V_R, f) \cdot f(V_R)
\end{aligned}
```

The indicator function $\phi$ can be defined in terms of a scaled activation function ($\sigma$, sigmoid) and the $P_v$ limits $(P_{vmin}, P_{vmax})$.
```math
\begin{aligned}
   \phi_L(V_R,f)&= \sigma(V_{rmin}-V_R)\sigma(-f)            \\
   \phi_U(V_R,f)&= \sigma(V_R-V_{rmax})\sigma(f)             \\
   \phi(V_R,f)  &= \left[1-\phi_L\right]\left[1-\phi_U\right]\\
\end{aligned}
```

The scale of the sigmoid function ($\alpha$ on the order of $10^3$) should be chosen so that for all practical parameters of the TGOV1 model, the sigmoid acts as a step function. This is further approximated by an algebraic form to obtain a practical function during implementation.
```math
\begin{aligned}
   \sigma(x) =\dfrac{1}{1+\exp(-\alpha x)}\approx \dfrac{1}{2}\left(\dfrac{\alpha x}{1+|\alpha x|}\right)
\end{aligned}
```

Lastly for the algebraic (continuous) piecewise functions, this implementation is straightforward when the approximation above is used.
```math
\begin{aligned}
    k_{sat}&=\sigma (E_{fd}' -S_A) \cdot S_B(E_{fd}' -S_A)^2   
\end{aligned}
```

## Initialization
At steady state we assume that $V_R$ is at or within its limits, 
and that the exciter is not saturated. The field $E_{fd}$ can be 
obtained through the steady-state conditions of the machine. 
We also assume for the moment that the stabilizer and over/under 
excitation limiters are non-existant.
```math
\begin{aligned}
    E_{fd}' &= E_{fd}  \\
    V_{f}  &= 0 \\
    V_{E}   &= 0 \\
    k_{sat} &= 0 \\
    V_{fx}   &= \dfrac{K_F}{T_F} E_{fd} \\
    V_R     &= K_E E_{fd} \\
    V_{tr}  &= \dfrac{1}{K_{a}}V_{R} \\
    E_C     &= V_{ref}-V_f - V_{tr}  \\
    V_{ts}  &= E_C \\
\end{aligned}
```
