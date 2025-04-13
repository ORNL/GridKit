# GENROU

## Block Diagram
<div align="center">
   <img align="center" src="../../../../../docs/Figures/GENROU.JPG">
   
   
  Figure 2: GENROU. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Simplifications
The GENROU model is a variation of the [General Synchronous Machine Model](../README.md)
- $`X''_{q}=X''_{d}`$
- $`X''_{d}`$ does not saturate
- Same relative amount of saturation occurs on both $`d`$ and $`q`$ axis


## Equations
### Nomenclature

See [Synchronous Machine Nomenclature](../README.md#nomenclature)

## Equations
### Differential Equations

See [Synchronous Machine Diffrential Equations](../README.md#differential-equations)

### Algebraic Equations

See [Synchronous Machine Algebraic Equations](../README.md#algebraic-equations)

## Initialization

### Without Saturation

Pressume there is no saturation to simplify solution procedure for initial conditions.

Using the power-flow solution, we have explicity solutions for the following variables. The internal variables $I_d$, $I_q$, $V_d$, and $V_q$ are calculated from the network interface equations. The remaining are algebraicillay solved from the steady-state initial conditions.

``` math
\begin{aligned}
\omega &= 0 \\
\delta &= \text{arg} \left[V_r + jV_i + (R_a + jX_q) (I_r + jI_i)\right] \\
  \psi^{''}_{d} &= V_q \\
  \psi^{''}_{q} &= -V_d \\
  \psi^{''} &= \sqrt{(\psi''_{d})^2+(\psi''_{q})^2} \\
  k_{sat}     &= S_B(\psi^{''}-S_A)^2 \\
  T_{elec}    &= (\psi''_{d} - I_dX_d^{''})I_q-(\psi''_{q} - I_qX_d^{''})I_d \\
  P_{mech} &= T_{elec} \\
\psi^{'}_d &=\dfrac{\psi^{''}_d/X_{d5}-X_{d2}I_d}{1+1/X_{d5}}\\
\psi^{'}_q &=\dfrac{X_{q2}I_q-\psi^{''}_q/X_{q5}}{1+1/X_{q5}}\\
E^{'}_d &=\psi^{'}_q - X_{q2}I_q \\
E^{'}_q &=\psi^{'}_d + X_{d2}I_d \\
E_{fd} &= E'_{q}+X_{d1}I_{d}+\psi^{''}_{d}k_{sat} \\
\end{aligned}
```

### With Saturation

It is important to point out that finding the initial value of $\delta$ for
the model without saturation direct method can be used. In case when saturation
is considered some "claver" math is needed. Key insight for determining initial
$\delta$ is that the magnitude of the saturation depends upon the magnitude
of $\psi''$, which is independent of $\delta$.

``` math
\begin{aligned}
  \delta=\tan^{-1}
  \left[
    \dfrac{(V_{i}+R_{a}I_{i})k_{sat}+(k_{sat}X''_{d}+X_{q}-X''_{q})I_{r}}
          {(V_{r}+R_{a}I_{r})k_{sat}-(k_{sat}X''_{d}+X_{q}-X''_{q})I_{i}}
  \right]
\end{aligned}
```
