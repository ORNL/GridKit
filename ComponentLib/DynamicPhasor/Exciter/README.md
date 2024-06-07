# **Exciter**


**Note: Exciter model not yet implemented**

<div align="center">
   <img align="center" src="/Documentation/Figures/EXDC1.JPG">
   
   
  Figure 1: Exciter EXDC1 model. Fifure courtesy of [PoweWorld](https://www.powerworld.com/WebHelp/).
</div>

## Nomenclature

### Inputs
- $`V_{REF}`$ - voltage reference set point
- $`E_{C}`$ - output from the terminal voltage transducer
- $`V_{S}`$ -  power system stabilizer output signal (if present)
- $`V_{UEL}`$ and $`V_{OEL}`$ - limiters

### States
- $`V_{t}`$ - terminal voltage (2 is sensed $`V_{t}`$)
- $`V_{B}`$ - input to a voltage regulator (3)
- $`V_{R}`$ - voltage regulator output also know as exciter field voltage (4)
- $`V_{F}`$ - stabilizing feedback signal (5)
### Parameters
- $`T_{R}`$ - filter time constant, sec (0)
- $`K_{A}`$ - voltage regulator gain (40)
- $`T_{A}`$ - time constant, sec (0.1)
- $`T_{B}`$ - lag time constant, sec (0)
- $`T_{C}`$ - lead time constant, sec (0)
- $`V_{RMAX}`$ - maximum control element output, pu (1)
- $`V_{RMIN}`$ - minimum control element output, pu (-1)
- $`K_{E}`$ - exciter field resistance line slope margine, pu (0.1)
- $`T_{E}`$ - exciter time constant, sec (0.5)
- $`K_{F}`$ - rate feedback gain, pu (0.05)
- $`T_{F1}`$ - rate feedback time constant, sec (0.7)
- $`E1`$ - field voltage value, 1 (2.8)
- $`SE1`$ - saturation factor at E1, (3.7)
- $`E2`$ - field voltage value, 2 (3.7)
- $`SE2`$ - saturation factor at E2, (0.33)

## Equations
First block
```math
\dfrac{dV_{t}}{dt}=\dfrac{1}{T_{R}}(E_{C}-V_{t})
```
Second block
```math
\dfrac{dx_{1}}{dt}=\dfrac{1}{T_{B}}((V_{REF}-V_{t}-V_{F}+V_{S}+V_{UEL}+V_{OEL})-V_{B})
```
```math
V_{B}=x_{1}+\dfrac{T_{C}}{T_{B}}(V_{REF}-V_{t}-V_{F}+V_{S}+V_{UEL}+V_{OEL})
```
Third block
```math
\dfrac{dV_{R}}{dt} = \begin{cases}
   \dfrac{1}{T_{A}}(K_{A}V_{B}-V_{R}) &\text{if } V_{RMIN}<=V_{R}<= V_{RMAX}\\
   0 &\text{if } V_{B}>0 \text{  and  } V_{R}>=V_{RMAX} &\text{ also then } V_{R}=V_{RMAX}\\
   0 &\text{if } V_{B}<0 \text{  and  } V_{R}<=V_{RMIN} &\text{ also then } V_{R}=V_{RMIN}\\
\end{cases}
```
Fourth block
```math
\dfrac{d\dfrac{E_{FD}}{\omega}}{dt}=\dfrac{1}{T_{E}}(V_{R}-\dfrac{(K_{E}+S_{E})E_{FD}}{\omega})
```
Feedback loop
```math
\dfrac{dx_{2}}{dt}=-\dfrac{V_{F}}{T_{F1}}
```
```math
V_{F}=x_{2}+\dfrac{K_{F}}{T_{F1}}\dfrac{E_{FD}}{\omega}
```
Saturation is modeled using an alternative quadratic function, with the value of Se specified at two points :
```math
Sat(x) = \begin{cases}
   \dfrac{B(x-A)^2}{x} &\text{if } x>A \\
   0 &\text{if } x<=A
\end{cases}
```
same as with the synchronous machines. There are two solutions, and one where $`A<1`$ should be chosen.
 
## Initialization
```math
V_{t}=V_{t_{0}}
```
```math
E_{C}=V_{t_{0}}
```
```math
(V_{REF}-V_{t}-V_{F}+V_{S}+V_{UEL}+V_{OEL})=V_{B}
```
```math
V_{R}=V{R_{0}}
```
```math
V_{B}=\dfrac{V{R_{0}}}{K_{A}}
```
```math
\dfrac{E_{FD}}{\omega}=\dfrac{E_{FD_{0}}}{\omega}
```
```math
V_{R}-\dfrac{(K_{E}+S_{E})E_{FD}}{\omega}=0
```
```math
V_{F}=0
```
```math
x_{2_{0}}=-\dfrac{K_{F}}{T_{F1}}\dfrac{E_{FD}}{\omega}


