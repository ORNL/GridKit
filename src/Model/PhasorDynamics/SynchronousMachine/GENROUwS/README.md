# GENROU

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
### State Variables

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
& k_{sat}       && \text{Saturation Coefficient} \\
& V_d, V_q      && \text{Machine Internal Voltage} \\
& T_{elec}      && \text{Electrical Torque} \\
& I_r, I_i      && \text{Terminal currents on the network reference frame} \\
& V_d, V_q      && \text{Terminal voltages on the network reference frame} \\
& I_d, I_q      && \text{Terminal currents on the machine d-q reference frame} \\
& P_{mech}      && \text{Mechanical power from the prime mover} \\
& E_{fd}        && \text{Field winding voltage from the excitation system} \\
\end{aligned}
```

## Equations



### Differential Equations

``` math
\begin{aligned}
  \dot\delta&=\omega\cdot\omega_0 \\
  \dot\omega&=\dfrac{1}{2H}\left(\dfrac{P_{mech}-D\omega}{1+\omega}-T_{elec}\right) \\


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
  \right) \\

  \dot{E}^{'}_{q} &= \dfrac{1}{T'_{d0}}
  \left(
    E_{fd}-E'_{q}-(X_{d}-X'_{d})
    \left(
      I_{d}+\dfrac{X'_{d}-X''_{d}}{(X'_{d}-X_{l})^2}(E'_{q}-\psi'_{d}-(X'_{d}-X_{l})I_{d})
    \right)
    -\psi''_{d}k_{sat}
  \right)\\
\end{aligned}
```

### Algebraic Equations
These algebraic equations define internal variables (7)
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

and the algebraic Network Interface Equations (4)
``` math
\begin{aligned}
\begin{cases}
  I_{dq}&=-jIe^{j\delta} \\
  V_{dq}&=-jVe^{j\delta}+ZI
\end{cases}
\implies
\begin{cases}
  I_d &= I_r \sin\delta - I_i\cos\delta \\
  I_q &= I_r \cos\delta + I_i\sin\delta \\

  V_d &= V_r \sin\delta - V_i\cos\delta + I_d R_a - I_q X^{''}_q\\
  V_q &= V_r \cos\delta + V_i\sin\delta + I_d X^{''}_q + I_q R_a\\
\end{cases} \\
\end{aligned}
```

## Simplifications

- $X^{''}_{q} = X^{''}_{d}$
- $X^{''}_{d}$ does not saturate
- The same relative amount of saturation occurs on both $d$ and $q$ axis



## Initialization

### Without Saturation

Pressume there is no saturation to simplify solution procedure for initial conditions. The variables () have explicity expressions, and the remainder are solved algebraically, shown below.

First, let the network-reference frame terminal variables be defined for notational purposes.

``` math
\begin{aligned}
\bar V &= V_r + jV_i && \text{Complex Terminal Voltage} \\
\bar I &= I_r + jI_i && \text{Complex Terminal Current} \\
\bar Z &= R_a + jX_q && \text{Machine Internal Impedance} \\
\end{aligned}
```

Using the power-flow solution, we have explicity solutions for the following variables. 

``` math
\begin{aligned}
\omega &= 0 \\
\delta &= \text{arg} (\bar V + \bar Z \bar I) \\
  I_{dq}&=-jIe^{j\delta} \\
  V_{dq}&=-jVe^{j\delta}+ZI \\
  \psi^{''}_{dq} &= jV_{dq} \\
  \psi^{''} &= |\psi^{''}_{dq}| \\
  k_{sat}     &= S_B(\psi^{''}-S_A)^2 \\
  T_{elec}    &= (\psi''_{d} - I_dX_d^{''})I_q-(\psi''_{q} - I_qX_d^{''})I_d \\
  P_{mech} &= T_{elec}
\end{aligned}
```

The remaining are algebraicillay solved from the steady-state initial conditions.

``` math
\begin{aligned}
\psi^{'}_d &= \dfrac{X^{'}_d-X_\ell}
  {1 + (X^{''}_d-X_\ell)}
  (\psi^{''}_d -(X^{''}_d-X_\ell)I_d ) \\
\psi^{'}_q &= \dfrac{X^{'}_q-X_\ell}
  {1 + (X^{''}_q-X_\ell)}
  (\psi^{''}_q -(X^{''}_q-X_\ell)I_q)\\
E^{'}_d &=\psi^{'}_q - (X^{'}_q-X_\ell)I_q \\
E^{'}_q &=\psi^{'}_d + (X^{'}_d-X_\ell)I_d \\
E_{fd} &= E'_{q}+(X_{d}-X'_{d})I_{d}+\psi^{''}_{d}k_{sat} \\
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
    \dfrac{(V_{iterm}+R_{a}I_{i})k_{sat}+(k_{sat}X''_{d}+X_{q}-X''_{q})I_{r}}
          {(V_{rterm}+R_{a}I_{r})k_{sat}-(k_{sat}X''_{d}+X_{q}-X''_{q})I_{i}}
  \right]
\end{aligned}
```

then

``` math
\begin{aligned}
    Sat(\psi'')=Sat(\vert V_{r}+jV_{i} \vert) 
\end{aligned}
```

where $k_{sat}$ is defined to be

$$
k_{sat}:=1+\left(\dfrac{X_{q}-X_{l}}{X_{d}-X_{l}}\right)Sat(\psi'')
$$

The following must be true (if not enforce the corrections):

``` math
\begin{aligned}
  X_{l}<=X''_{q}<=X'_{q}<=X_q \\
  X_{l}<=X''_{d}<=X'_{d}<=X_d
\end{aligned}
```
