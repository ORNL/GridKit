# GENSAL 
## Simplifications

- $`X''_{q}=X''_{d}`$
- $`X''_{d}`$ does not saturate
- only $`d`$ axis affected by saturation
- $`X_{q}=X'_{q}`$
- $`T'_{q0}`$ is neglected

<div align="center">
   <img align="center" src="/Documentation/Figures/GENSAL.JPG">
   
   
  Figure 2: GENSAL. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Equations
### Algebraic Equations


- Fluxes

``` math
E''_{d}=\psi''_{q}
```
``` math
E''_{q}=\psi''_{d}=+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}}
```
```math
\psi_{d}=-I_{d}X''_{d}+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}}=-I_{d}X''_{d}+E''_{q}
```
```math
\psi_{q}=-I_{q}X''_{d}+\psi''_{q}
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
T'_{d0}\dfrac{dE'_{q}}{dt}=E_{fd}-E'_{q}-(X_{d}-X'_{d})(I_{d}-\dfrac{X'_{d}-X''_{d}}{(X'_{d}-X_{l})^2}(+\psi'_{d}+(X'_{d}-X_{l})I_{d}-E'_{q}))-\psi''_{d}Sat(\psi'')
```
```math
T''_{d0}\dfrac{d\psi'_{d}}{dt}=-\psi'_{d}-(X'_{d}-X_{l})I_{d}+E'_{q}
```
```math
T''_{q0}\dfrac{d\psi''_{q}}{dt}=-\psi''_{q}-(X_{q}-X''_{q})I_{q}
```






























