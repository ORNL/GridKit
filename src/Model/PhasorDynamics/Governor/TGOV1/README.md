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
- $V_{max}$ - Max Valve Position
- $V_{min}$ -  Min Valve Position
- $T_2$ - Turbine Numerator Time Constat
- $T_3$ - Turbine Delay
- $D_t$ - Turbine Damping Coefficient

## Equations


### Algebraic Equations

```math
\Delta\omega=\dfrac{\omega-\omega_{s}}{\omega_{s}}
```

### Differential Equations
First block
```math
\dfrac{dV}{dt} = \begin{cases}
   \dfrac{1}{T1}(\dfrac{P_{REF}-\Delta\omega}{R}-V) &\text{if } V_{min}<=V<= V_{max}\\
   0 &\text{if } \dfrac{P_{REF}-\Delta\omega}{R}>0 \text{  and  } V>=V_{max} &\text{ also then } V=V_{max}\\
   0 &\text{if } \dfrac{P_{REF}-\Delta\omega}{R}<0 \text{  and  } V<=V_{min} &\text{ also then } V=V_{min}\\
\end{cases}
```
Second block
```math
\dfrac{dx}{dt}=\dfrac{1}{T3}(V-P)
```
```math
P=x+\dfrac{T2}{T3}V
```
Output
```math
P_{mech}=P-\Delta\omega D_{t}
```
