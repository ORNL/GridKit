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
### Auxillary Parameters
``` math
\begin{aligned}
X_{d1} &= X_d-X_d{'}      & X_{q1} &= X_q-X_q{'} \\
X_{d2} &= X_d^{'}-X_\ell  & X_{q2} &= X_q^{'}-X_\ell\\
X_{d3} &= (X_d^{'}-X_d^{''})/X_{d2}^2 & X_{q3} &= (X_q^{'}-X_q^{''})/X_{q2}^2 \\
X_{d4} &= (X_d^{'}-X_d^{''})/X_{d2}   & X_{q4} &= (X_q^{'}-X_q^{''})/X_{q2} \\
X_{d5} &= (X_d^{''}-X_\ell)/X_{d2}    & X_{q5} &= (X_q^{''}-X_\ell)/X_{q2}\\
X_{qd} &= (X_q-X_\ell)/(X_d-X_\ell)
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
  \dot\omega&=\dfrac{1}{2H}\left(\dfrac{P_{mech}-D\omega}{1+\omega}-T_{elec}\right)\\
  \dot{\psi}^{'}_{d} &= \dfrac{1}{T''_{d0}}(E'_{q}-\psi'_{d}-X_{d2}I_{d})\\
  \dot{\psi}^{'}_{q} &= \dfrac{1}{T''_{q0}}(E'_{d}-\psi'_{q}+X_{q2}I_{q})\\


  \dot{E}^{'}_{d} &= \dfrac{1}{T'_{q0}}
  \left(  -E'_{d}+X_{q1}
    (I_{q}-X_{q3}(E'_{d}-\psi'_{q}+X_{q2}I_{q}))
    + X_{qd}\psi''_{q}k_{sat}
  \right) \\

  \dot{E}^{'}_{q} &= \dfrac{1}{T'_{d0}}
  \left(
    E_{fd}-E'_{q}-X_{d1}
    (I_{d}+X_{d3}(E'_{q}-\psi'_{d}-X_{d2}I_{d}))
    -\psi''_{d}k_{sat}
  \right)\\
\end{aligned}
```

### Algebraic Equations

These algebraic equations define internal variables (7)
``` math
\begin{aligned}
  \psi''_{q} &= -E'_{d}X_{q5}
                -\psi'_{q}X_{q4} \\
  \psi''_{d}  &= +E'_{q}X_{d5}
                +\psi'_{d}X_{d4}\\
  \psi^{''}   &= \sqrt{(\psi''_{d})^2+(\psi''_{q})^2} \\
  k_{sat}     &= S_B(\psi^{''}-S_A)^2 \\
  V_{d}       &= -\psi''_{q}(1+\omega)\\
  V_{q}       &= +\psi''_{d}(1+\omega)\\
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
  \psi^{''} &= |\psi^{''}_{dq}| \\
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
    \dfrac{(V_{iterm}+R_{a}I_{i})k_{sat}+(k_{sat}X''_{d}+X_{q}-X''_{q})I_{r}}
          {(V_{rterm}+R_{a}I_{r})k_{sat}-(k_{sat}X''_{d}+X_{q}-X''_{q})I_{i}}
  \right]
\end{aligned}
```
