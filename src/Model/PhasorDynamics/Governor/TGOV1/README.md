# **TGOV1**

## Control Diagram

Standard model of the stream turbine

<div align="center">
   <img align="center" src="../../../../../docs/Figures/TGOV1.JPG">
   
   
  Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Nomenclature

### Algebraic Variables
- $\omega$ - Per-unit Speed Deviation
- $P_{ref}$ - Reference Power
- $P_{mech}$ - Mechnical Power to Generator

### Differential Variables
- $P_{tx}$ - Turbine Power 
- $P_{v}$ - Valve Position

### Parameters
- $R$ - Droop Constant
- $T_1$ - Valve Time Delay
- $T_2$ - Turbine Numerator Time Constat
- $T_3$ - Turbine Delay
- $V_{max}$ - Max Valve Position
- $V_{min}$ -  Min Valve Position
- $D_t$ - Turbine Damping Coefficient

## Equations


### Algebraic Equations
The algebraic equation dictating the mechnical power output.
```math
\begin{aligned}
   P_{mech} &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \omega \\
\end{aligned}
```

### Differential Equations
The TGOV1 differential equations, which include a piecewise definition.
The domain of the state variable $P_{v}\in[V_{min}, V_{max}]$ is enforced
through the piece-wise definition below.
```math
\begin{aligned}
   \dot{P}_{tx}   &= P_v - \dfrac{1}{T_3}(P_{tx}+T_2P_v) \\
   \dot{P}_{v}    &= 
      \dfrac{1}{T_1}
   \begin{cases}
      -P_{v} + \frac{1}{R}(P_{ref}-\omega)
         &  \text{if } P_v \in \in[V_{min}, V_{max}]\\
      0  
         &  \text{else } \\
   \end{cases}
\end{aligned}
```
