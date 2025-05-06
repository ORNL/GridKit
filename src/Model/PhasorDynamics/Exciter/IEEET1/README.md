# **IEEET1**

## Control Diagram

Standard model of the IEEET1 Exciter.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_IEEET1_Diagram.png">
   
   
  Figure 1: Exciter IEEET1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Nomenclature

### Algebraic Variables
- $V_{tr}$ - Terminal Voltage Error input to controller
- $V_{F}$ - Feedback Voltage
- $V_{E}$ - Excitation control voltage
- $E_{fd}$ - Field winding voltage
- $k_{sat}$ - Saturation variable

These can be constants or external states
- $E_{C}$ - Compensated machine terminal voltage magnitude
- $V_{ref}$ - Referebce ternubak voltage
- $V_{UEL}$ - Input from under excitation limiter
- $V_{OEL}$ - Input from over excitation limiter
- $V_{S}$ - Input from stabilizer controller
- $\omega$ - Machine speed deviation from machine model

### Differential Variables
- $V_{ts}$ - Sensed terminal voltage
- $V_{R}$ - Voltage regulator
- $E_{fd}'$ - Field-current pre-speed multiplier
- $V_{fx}$ - Exciter feedback internal state

### Parameters
- $T_R$ - Time constant for voltage sensing
- $K_a, T_a$ - Coefficient and time constant for voltage regulation
- $K_e, T_e$ - Coefficient and time constant for excitation system
- $K_f, T_f$ - Coefficient and time constant for feedback
- $V_{rmin}, V_{rmax}$ - Limits to voltage regulation
- $E_1, S_{e1}, E_{2}, S_{e2}$ - Saturation Parameters
- $I_{spdlm}$ - Speed Limit flag indicator

## Equations


### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   V_{tr} &= V_{ref} - V_{ts}+V_{UEL} + V_{OEL} + V_S - V_F\\
    V_{f} &= \dfrac{E_{fd}' K_F}{T_F} - V_{fx}\\
    E_{fd}&= \begin{cases}
        E_{fd}^{'}           &  \text{if } I_{spdlm}\\
        (1+\omega)E_{fd}^{'}  &  \text{else } \\
   \end{cases}\\
    k_{sat}&= \begin{cases}
        S_B(E_{fd}^{'} -S_A)^2        &  \text{if } E_{fd}^{'} >S_A\\
        0  &  \text{else } \\
   \end{cases} \\
    V_{E} &= k_{sat}\cdot E_{fd}^{'} \\
\end{aligned}
```


### Differential Equations
The TGOV1 differential equations, as derived from the model diagram.
```math
\begin{aligned}
   \dot{V}_{ts}   &= \dfrac{1}{T_R}(E_C-V_{ts}) \\
   \dot{V}_{R}    &= 
      \dfrac{1}{T_A}
   \begin{cases}
      -V_{R}+K_{a}V_{tr}
         &  \text{if } V_R \in (V_{rmin}, V_{rmax})\\
      0  
         &  \text{else } \\
   \end{cases}
\end{aligned}
```
The domain of the state variable $V_{R}\in(V_{rmin}, V_{rmax})$ is enforced
through the piece-wise definition above. However, depending on the
general solver's requirments, this may need to be expressed as a
smooth approximation (bump function/smooth indicator). Perhaps like the 
following.
```math
\begin{aligned}
   I(x,x_0,\epsilon) 
            &= \dfrac{1}{2} + \dfrac{1}{2}\tanh 
            \left(\dfrac{x-x_0}{\epsilon}\right) \\
   \phi(x)  &= I \left(x;V_{rmin}, \epsilon \right) -
            I \left(x;V_{rmax}, \epsilon \right)
            \qquad \epsilon <<1 \\
   \dot{V}_{R}   
            &= 
            \dfrac{1}{T_A}\cdot \phi(V_R)
            \left[
               -V_{R}+K_{a}V_{tr}
            \right] \\
\end{aligned}
```
