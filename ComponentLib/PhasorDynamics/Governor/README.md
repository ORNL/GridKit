# **Governor**

## TGOV1 Model


**Note: Governor model not yet implemented**

Standard model of the stream turbine

<div align="center">
   <img align="center" src="/Documentation/Figures/TGOV1.JPG">
   
   
  Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Nomenclature

### Inputs 
- $`P_{REF}`$ - reference power (set point)
- $`\Delta\omega`$
### States 
- $`P`$ - turbine power (state 1 in the Figure)
- $`V`$ - valve position (state 2 in the Figure)
### Parameters
- $`R`$ - permanent droop, pu (0.05)
- $`T2`$ - steam bowl time constant, sec (0.5)
- $`V_{max}`$ - maximum valve position limit (1)
- $`V_{min}`$ -  minimum valve position limit (0)
- $`T2`$ - numerator time constant of T2, T3 block, sec (2.5)
- $`T3`$ - reheater time constant, sec (7.5)
- $`D_{t}`$ - turbine damping coefficient, pu (0)
- $`T_{rate}`$ - turbine rating, MW (0) 

## Equations
```math
\Delta\omega=\dfrac{\omega-\omega_{s}}{\omega_{s}}
```
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