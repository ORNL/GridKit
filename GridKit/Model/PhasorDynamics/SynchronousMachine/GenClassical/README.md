# Classical Generator

An electrical machine model with two differential variables (i.e. second-order
model) is often called classical generator model. While its predictive ability
is limited, it is useful for studies of grid network properties. Mathematically,
it is equivalent to a driven damped pendulum model.

## Model Parameters

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------|----------------------
$\omega_0$  | [rad/s] | synchronous frequency           |
$H$         | [s]     | rotor inertia                   |
$D$         | [p.u.]  | damping coefficient             |
$R_a$       | [p.u.]  | winding resistance              |
$X_{dp}$    | [p.u.]  | machine reactance parameter     |

### Model Derived Parameters

- $G = \dfrac{R_a}{R_a^2 + X_{dp}^2} ~~~$ equivalent stator winding conductance
- $B = \dfrac{-X_{dp}}{R_a^2 + X_{dp}^2} ~~~$ equivalent stator winding susceptance

<br>

## Model Variables

### Internal Variables

#### Differential

Symbol      | Units   | Description         | Note
------------|---------|---------------------|----------------------
$\delta$    | [rad]   | machine power angle |
$\omega$    | [p.u]   | machine speed deviation       | Optionally read by a governor or a stabilizer component

#### Algebraic

Symbol  | Units  | Description                         | Note
--------|--------|-------------------------------------|-------------
$T_{e}$ | [p.u.] | electrical torque                   |
$I_r$   | [p.u.] | machine real injection current      | read by bus
$I_i$   | [p.u.] | machine imaginary injection current | read by bus

Note: All three can be expressed as a function called by the model equations. We add
these as variables as they are needed for outputs.

<br>

### External Variables

External variables enter component model equations but are owned by other
components. The other components also provide equations needed to have a
balanced system of equations. 

#### Differential

None.

#### Algebraic

Symbol | Units   | Description                   | Note
-------|---------|-------------------------------|----------------------
$V_r$  | [p.u.]  | machine bus real voltage      | owned by a bus object
$V_i$  | [p.u.]  | machine bus imaginary voltage | owned by a bus object
$P_m$  | [p.u.]  | mechanical power input        | owned by governor, constant if no governor is connected to the machine
$E_p$  | [p.u.]  | field winding voltage         | owned by exciter, constant if no exciter is connected to the machine

<br>


## Model Equations

### Differential Equations

```math 
\begin{aligned}
\dot{\delta} &= \omega \cdot \omega_0 \\
\dot{\omega} &= \frac{1}{2H}\left( \frac{P_{m} - D\omega}{1+\omega}   - T_{e}\right)
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
    0 &= T_{e} - \left( G E_p^2 - E_p \left[(G V_r - B V_i)\cos\delta + (B V_r + G V_i)\sin\delta \right]\right) \\
    0 &= I_r + G V_r - B V_i - E_p(G \cos\delta - B \sin\delta) \\
    0 &= I_i + B V_r + G V_i - E_p(B \cos\delta + G \sin\delta)
\end{aligned}
```
As noted earlier, all three algebraic equations can be expressed as functions
and substituted directly in the component and bus equations, respectively. We
use redundant variables for modeling convenience.

<br>

## Initialization

To initialize the model, given bus voltages $V_r$, $V_i$, and initial generator
injection active and reactive power, $P$ and $Q$, we take following steps to
initialize the system:

Complex power is defined as
```math
S=VI^{*}
```
or
```math
P + jQ = (V_r + j V_i)(I_r - j I_i).
```
From here, we compute injection currents from the initial power injection and bus
voltages as
```math
\begin{aligned}
I_r &= \frac{PV_r + QV_i}{V_r^2 + V_i^2} \\
I_i &= \frac{PV_i - QV_r}{V_r^2 + V_i^2}
\end{aligned} 
```

We substitute the expressions above into equations for current injections and
obtain
```math
\begin{aligned}
E_p \sin\delta &= \dfrac{-B I_r + G I_i}{G^2 + B^2} + V_i \\
E_p \cos\delta &= \dfrac{G I_r + B I_i}{G^2 + B^2} + V_r
\end{aligned}
```
By dividing these two equations, we get an expression for the machine angle at the
steady state:
```math
\delta = \arctan \dfrac{E_i}{E_r} \, ,
```
And by squaring and adding them, we get an expression for the field
winding voltage at the steady state
```math
E_p = \sqrt{E_r^2 + E_i^2}   \, ,
```
where
```math
\begin{aligned}
E_r &= \dfrac{G I_r + B I_i}{G^2 + B^2} + V_r \, ,\\
E_i &= \dfrac{-B I_r + G I_i}{G^2 + B^2} + V_i \, .
\end{aligned}
```

Next, we set the machine speed deviation to zero:
```math
\omega = 0
```

Now, we can compute the electrical torque and set the mechanical torque to be equal
to the electrical:
```math
\begin{aligned}
T_{e} &= G E_p^2 - E_p \left[ (G V_r - B V_i ) \cos\delta + (B V_r + G V_i )\sin\delta \right] \\
P_{m} &= T_{e} 
\end{aligned} 
```

With this, we initialize the machine at a steady state.


## Model Outputs

Real and imaginary currents, $I_r$ and $I_i$, are model algebraic variables,
oriented as leaving the machine (i.e. entering the bus).

Machine angle $\delta$ and speed $\omega$ are model state variables.

Active and reactive power ($P$, and $Q$) 
are the real and imaginary parts of the complex power, 
defined as $S=VI^{\ast}=(V_r + j V_i)(I_r - jI_i)$
``` math
\begin{aligned}
      P &= V_{r} I_{r} + V_{i} I_{i}\\
      Q &= V_{i} I_{r} - V_{r} I_{i}
\end{aligned}
```
Power outputs are oriented leaving the machine (i.e. entering the bus).