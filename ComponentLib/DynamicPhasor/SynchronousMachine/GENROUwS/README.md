# GENROU
## Simplifications

- $`X''_{q}=X''_{d}`$
- $`X''_{d}`$ does not saturate
- same relative amount of saturation occurs on both $`d`$ and $`q`$ axis

<div align="center">
   <img align="center" src="/Documentation/Figures/GENROU.JPG">
   
   
  Figure 2: GENROU. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

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
T'_{d0}\dfrac{dE'_{q}}{dt}=E_{fd}-E'_{q}-(X_{d}-X'_{d})(I_{d}-\dfrac{X'_{d}-X''_{d}}{(X'_{d}-X_{l})^2}(+\psi'_{d}+(X'_{d}-X_{l})I_{d}-E'_{q}))-\psi''_{d}Sat(\psi'')
```
```math
T''_{d0}\dfrac{d\psi'_{d}}{dt}=-\psi'_{d}-(X'_{d}-X_{l})I_{d}+E'_{q}
```
```math
T''_{q0}\dfrac{d\psi'_{q}}{dt}=-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d}
```
```math
T'_{q0}\dfrac{dE'_{d}}{dt}= -E'_{d}+(X_{q}-X'_{q})(I_{q}-\dfrac{X'_{q}-X''_{q}}{(X'_{q}-X_{l})^2}(-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d}))+\psi''_{q}(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}})Sat(\psi'')
```
## Initialization

From the block diagram it can be written:

```math
-\psi'_{d}-(X'_{d}-X_{l})I_{d}+E'_{q}=0
```
``` math
-\psi''_{d}+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}}=0
```
```math
-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d}=0
```
``` math
\psi''_{q}+E'_{d}\dfrac{X''_{q}-X_{l}}{X'_{q}-X_{l}}+\psi'_{q}\dfrac{X'_{q}-X''_{q}}{X'_{q}-X_{l}}=0
```
```math
-E'_{d}+(X_{q}-X'_{q})I_{q}+\psi''_{q}(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}})Sat(\psi'')=0
```

Internal voltage on the referece frame can be calculated directly:
```math
V_{r}=V_{rterm}+R_{a}I_{r}-X''_{d}I_{i}
```
``` math
V_{i}=V_{iterm}+R_{a}I_{i}-X''_{d}I_{r}
```
then
```math
Sat(\psi'')=Sat(\vert V_{r}+jV_{i} \vert)
```

It is important to point out that finding the initial value of $`\delta`$ for the model without saturation direct method can be used. In case when saturation is considered some "claver" math is needed. Key insight for determining initial $`\delta`$ is that the magnitude of the saturation depends upon the magnitude of $`\psi''`$, which is independent of $`\delta`$.
```math
\delta=tan^{-1}(\dfrac{K_{sat}V_{iterm}+K_{sat}R_{a}I_{i}+(K_{sat}X''_{d}+X_{q}-X''_{q})I_{r}}{K_{sat}V_{rterm}+K_{sat}R_{a}I_{r}-(K_{sat}X''_{d}+X_{q}-X''_{q})I_{i}})
```
where
```math
K_{sat}=(1+(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}})Sat(\psi''))
```
Following must be true (if not enforce the corrections):

```math
X_{l}<=X''{q}<=X'{q}<=Xq
```
```math
X_{l}<=X''{d}<=X'{d}<=Xd
```
