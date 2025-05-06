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

```math
\begin{aligned}
   \dot{P}_{mech} &= \dfrac{1}{T_3}(P_{tx}+T_2P_v) - D_t \omega \\
\end{aligned}
```

### Differential Equations
First block
```math
\begin{aligned}
   \dot{P}_{tx}   &= P_v - \dfrac{1}{T_3}(P_{tx}+T_2P_v) \\
   \dfrac{dV}{dt} &= 
   \begin{cases}
      \dfrac{1}{T1}(\dfrac{P_{ref}-\Delta\omega}{R}-V) 
         &  \text{if } P_v \geq V_{max}\\
      0  &  \text{if } P_v \leq V_{min}\\
         &  \text{else } \\
   \end{cases}
\end{aligned}
```
