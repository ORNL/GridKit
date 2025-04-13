# GENROU
## Simplifications

- $`X''_{q}=X''_{d}`$
- $`X''_{d}`$ does not saturate
- same relative amount of saturation occurs on both $`d`$ and $`q`$ axis

<div align="center">
   <img align="center" src="../../../../../docs/Figures/GENROU.JPG">
   
   
  Figure 2: GENROU. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Nomenclature
### Parameters
``` math
\begin{aligned}
& \omega_0  && \text{Nominal Frequnecy}\\
& H         && \text{Inertia Constant} \\
& D         && \text{Damping Coefficient}\\
& R_a       && \text{Winding Resistance}\\
& T^{'}_{d0}, T^{''}_{d0}, T^{'}_{q0}, T^{''}_{q0} 
            && \text{Time Constants} \\
& X_d, X^{'}_d, X^{''}_d,  X_q, X^{'}_q, X^{''}_q, X_\ell 
            && \text{Machine Reactance Parameters} \\
& S_{10}, S_{12} \to S_{A}, S_{B} && \text{Saturation Parameters}
\end{aligned}
```

### States Variables

``` math
\begin{aligned}
& \delta        && \text{Machine Internal Angle}\\
& \omega        && \text{Machine Relative Speed} \\
& E^{'}_q,\psi^{'}_d, \psi^{'}_q, E^{'}_d         
                && \text{Machine Internal Flux Values}\\
& \psi^{''}_q, \psi^{''}_d, \psi^{''}       
                && \text{Machine Total Subtransient Flux}\\
\end{aligned}
```

### Algebraic Variables

``` math
\begin{aligned}
& k_{sat}             
                && \text{Saturation Coefficient} \\
& V_d, V_q        
                && \text{Machine Internal Voltage} \\
& T_{elec} 
                && \text{Electrical Torque} \\
& I_r, I_i 
                && \text{Terminal currents on the network reference frame} \\
& V_d, V_q 
                && \text{Terminal voltages on the network reference frame} \\
& I_d, I_q 
                && \text{Terminal currents on the machine d-q reference frame} \\
& P_{mech} 
                && \text{Mechanical power from the prime mover} \\
& E_{fd} 
                && \text{Field winding voltage from the excitation system} \\
\end{aligned}
```

## Equations



### Differential Equations

``` math
\begin{aligned}
  \dot\delta&=\omega\cdot\omega_0 \\
  \dot\omega&=\dfrac{1}{2H}\left(\dfrac{P_{mech}-D\omega}{1+\omega}-T_{elec}\right) \\
  \dot{E}^{'}_{q} &= \dfrac{1}{T'_{d0}}
  \left(
    E_{fd}-E'_{q}-(X_{d}-X'_{d})
    \left(
      I_{d}+\dfrac{X'_{d}-X''_{d}}{(X'_{d}-X_{l})^2}(E'_{q}-\psi'_{d}-(X'_{d}-X_{l})I_{d})
    \right)
    -\psi''_{d}k_{sat}
  \right)\\

  \dot{\psi}^{'}_{d} &= \dfrac{1}{T''_{d0}}
  \left(
    E'_{q}-\psi'_{d}-(X'_{d}-X_{l})I_{d}
  \right)\\

  \dot{\psi}^{'}_{q} &= \dfrac{1}{T''_{q0}}
  \left(
    E'_{d}-\psi'_{q}+(X'_{q}-X_{l})I_{q}
  \right)
  \\


  \dot{E}^{'}_{d} &= \dfrac{1}{T'_{q0}}
  \left(  -E'_{d}+(X_{q}-X'_{q})
    \left(
      I_{q}-\dfrac{X'_{q}-X''_{q}}{(X'_{q}-X_{l})^2}(E'_{d}-\psi'_{q}+(X'_{q}-X_{l})I_{q})
    \right)
    + 
    \left( 
      \dfrac{X_{q}-X_{l}}{X_{d}-X_{l}} 
    \right)
  \psi''_{q}k_{sat}
  \right)
\end{aligned}
```

### Algebraic Equations
These algebraic equations define internal variables
``` math
\begin{aligned}
  -\psi''_{q} &= E'_{d}\dfrac{X''_{q}-X_{l}}{X'_{q}-X_{l}}
                +\psi'_{q}\dfrac{X'_{q}-X''_{q}}{X'_{q}-X_{l}} \\
  \psi''_{d}  &= E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}
                +\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}} \\
  \psi^{''}   &= \sqrt{(\psi''_{d})^2+(\psi''_{q})^2} \\
  k_{sat}     &= S_B(\psi^{''}-S_A)^2 \\
  V_{d}       &= -\psi''_{q}(1+\omega)\\
  V_{q}       &= \psi''_{d}(1+\omega)\\
  T_{elec}    &= (\psi''_{d} - I_dX_d^{''})I_q-(\psi''_{q} - I_qX_d^{''})I_d
\end{aligned}
```

and the algebraic Network Interface Equations
``` math
\begin{aligned}
  I_d &= I_r \sin\delta - I_i\cos\delta \\
  I_q &= I_r \cos\delta + I_i\sin\delta \\
  V_d &= V_r \sin\delta - V_i\cos\delta + I_d R_a - I_q X^{''}_q\\
  V_q &= V_r \cos\delta + V_i\sin\delta + I_d X^{''}_q + I_q R_a\\
\end{aligned}
```

## Initialization

From the block diagram it can be written:

``` math
\begin{aligned}
  0&=-\psi'_{d}-(X'_{d}-X_{l})I_{d}+E'_{q} \\
  0&=-\psi''_{d}+E'_{q}\dfrac{X''_{d}-X_{l}}{X'_{d}-X_{l}}+\psi'_{d}\dfrac{X'_{d}-X''_{d}}{X'_{d}-X_{l}} \\
  0&=-\psi'_{q}+(X'_{q}-X_{l})I_{q}+E'_{d} \\
  0&=\psi''_{q}+E'_{d}\dfrac{X''_{q}-X_{l}}{X'_{q}-X_{l}}+\psi'_{q}\dfrac{X'_{q}-X''_{q}}{X'_{q}-X_{l}} \\
  0&= -E'_{d}+(X_{q}-X'_{q})I_{q}+(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}})\psi''_{q}k_{sat}
\end{aligned}
```

Internal voltage on the referece frame can be calculated directly:
``` math
\begin{aligned}
V_{r}&=V_{rterm}+R_{a}I_{r}-X''_{d}I_{i} \\
V_{i}&=V_{iterm}+R_{a}I_{i}-X''_{d}I_{r}
\end{aligned}
```

then

``` math
\begin{aligned}
    Sat(\psi'')=Sat(\vert V_{r}+jV_{i} \vert) 
\end{aligned}
```

It is important to point out that finding the initial value of $`\delta`$ for
the model without saturation direct method can be used. In case when saturation
is considered some "claver" math is needed. Key insight for determining initial
$`\delta`$ is that the magnitude of the saturation depends upon the magnitude
of $`\psi''`$, which is independent of $`\delta`$.

``` math
\begin{aligned}
  \delta=\tan^{-1}
  \left(
    \dfrac{K_{sat}V_{iterm}+K_{sat}R_{a}I_{i}+(K_{sat}X''_{d}+X_{q}-X''_{q})I_{r}}
          {K_{sat}V_{rterm}+K_{sat}R_{a}I_{r}-(K_{sat}X''_{d}+X_{q}-X''_{q})I_{i}}
  \right)
\end{aligned}
```

where

$$
K_{sat}=(1+(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}})Sat(\psi''))
$$

Following must be true (if not enforce the corrections):

``` math
\begin{aligned}
  X_{l}<=X''{q}<=X'{q}<=Xq \\
  X_{l}<=X''{d}<=X'{d}<=Xd
\end{aligned}
```
