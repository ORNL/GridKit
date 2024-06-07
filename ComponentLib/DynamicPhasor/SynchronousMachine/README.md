# **Synchronous Machine - GENERAL**


**Note: Synchronous machine models not yet implemented**



## Convention


<div align="center">
   <img align="center" src="/Documentation/Figures/SM1.JPG">
   
   
  Figure 1: Synchronous Machine. Figure courtesy of [PowerWorld](https://www.powerworld.com/files/Synchronous-Machines.pdf)
</div>

q-axis leads the d-axis

rotor angle w.r.t. to q-axis

## Types

Two main types:

- Round Rotor (we will use GENROU model)
- Salient Rotor/ Salient Pole (we will use GENSAL model)

## Nomenclature
### Variables
- $`\delta`$ - rotor angle
- $`\omega_{s}`$ - synchronous speed (2$`\pi`$60)
- $`\Delta \omega_{pu}`$ - deviation of rotor speed away from synchronous speed
- $`I_{d}, I{q},I_{0}`$ - stator currents
- $`V_{dterm}, V_{qterm}, V_{0term}`$ - stator voltages
- $`\psi_{d}, \psi_{q}, \psi_{0}`$ - stator flux
- $`E'_{d}, E'_{q}, \psi'_{q}, \psi'_{d}`$ - rotor fluxes
- $`E_{fd}`$ - field voltage input (from exciter)
- $`P_{mech}`$ - mechanical
### Parameters
- $`H`$ - intertia constant, sec (3)
- $`D`$ - damping factor, pu (0)
- $`R_{s}`$ - stator resistance, pu (0)
- $`X_{l}`$ - stator leakage reactance, pu (0.15)
- $`X_{d}`$ - direct axis synchronous reactance, (2.1)
- $`X'_{d}`$ - direct axis transient reactance, (0.2)
- $`X''_{d}`$ - direct axis sub-transient reactance, (0.18)
- $`X_{q}`$ - quadrature axis synchronous reactance, (0.5)
- $`X'_{q}`$ - quadrature axis transient reactance, (0.47619)
- $`X''_{q}`$ - quadrature axis sub-transient reactance, (0.18)
- $`T'_{d0}`$ - open circuit direct axis transient time const., (7)
- $`T''_{d0}`$ - open circuit direct axis sub-transient time const., (0.04)
- $`T'_{q0}`$ - open circuit quadrature axis transient time const., (0.75)
- $`T''_{q0}`$ - open circuit quadrature axis sub-transient time const., (0.05) 
- $`S1`$ - saturation factor at 1.0 pu flux, (0) 
- $`S12`$ - saturation factor at 1.2 pu flux, (0) 

## Equations

### Algebraic Equations


- Fluxes

``` math
E''_{d}=-\psi''_{q}=+E'_{d}\dfrac{X''_{q}-X_{l}}{X'_{q}-X_{l}}+\psi'_{q}\dfrac{X'_{q}-X''_{q}}{X'_{q}-X_{l}}
```
``` math
E''_{q}=\psi''_{d}=+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}}
```
```math
\psi_{d}=-I_{d}X''_{d}+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}}=-I_{d}X''_{d}+E''_{q}
```
```math
\psi_{q}=-I_{q}X''_{q}-E'_{d}\dfrac{X''_{q}-X_{l}}{X'_{q}-X_{l}}-\psi'_{q}\dfrac{X'_{q}-X''_{q}}{X'_{q}-X_{l}}=-I_{q}X''_{q}-E''_{d}
```
- Stator
``` math
V_{dterm}=E''_{d}(1+\Delta\omega_{pu})-R_{s}I_{d}+X''_{q}I_{q}
```
``` math
V_{qterm}=E''_{q}(1+\Delta\omega_{pu})-R_{s}I_{q}-X''_{d}I_{d}
```

### Differential Equations


- Mechanical Dynamic Equations
``` math
\dfrac{d\delta}{dt}=\Delta \omega_{pu}*\omega_{s}
```
``` math
2H\dfrac{d\omega}{dt}=\dfrac{P_{mech}-D\omega}{1+\Delta\omega_{pu}}-(\psi_{d}I_{q}-\psi_{q}I_{d})
```
- Rotor Dynamic Equations
```math
T'_{d0}\dfrac{dE'_{q}}{dt}=E_{fd}-E'_{q}-(X_{d}-X'_{d})(I_{d}-\dfrac{X'_{d}-X''_{d}}{(X'_{d}-X_{l})^2}(+\psi'_{d}+(X'_{d}-X_{l})I_{d}-E'_{q}))
```
```math
T''_{d0}\dfrac{d\psi'_{d}}{dt}=-\psi'_{d}-(X'_{d}-X_{l})I_{d}+E'_{q}
```
```math
T''_{q0}\dfrac{d\psi'_{q}}{dt}=-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d}
```
```math
T'_{q0}\dfrac{dE'_{d}}{dt}= -E'_{d}+(X_{q}-X'_{q})(I_{q}-\dfrac{X'_{q}-X''_{q}}{(X'_{q}-X_{l})^2}(-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d}))
```
Previos equations can be used to model any machine, however ***SATURATION*** is missing.

Saturation means increasingly large amounts of current are needed to increase the flux density. There are various methods to include the saturation (it is not standardized yet). We are going to use the approach implemented in PTI PSSS/E and PowerWorld Simulator (scaled quadratic). 

```math
Sat(x) = \begin{cases}
   \dfrac{B(x-A)^2}{x} &\text{if } x>A \\
   0 &\text{if } x<=A
\end{cases}
```
There are two solutions, and one where $`A<1`$ should be chosen.

Hint!

Negative values are not allowed.
