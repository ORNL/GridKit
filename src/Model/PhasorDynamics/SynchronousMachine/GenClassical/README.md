# Classical Generator

An electrical machine model with two differential variables (i.e. second order
model) is often called classical generator model. While its predicitve ability
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

- $g = \dfrac{R_a}{R_a^2 + X_{dp}^2}$
- $b = \dfrac{-X_{dp}}{R_a^2 + X_{dp}^2}$

<br>

## Model Variables

### Internal Variables

#### Differential

Symbol      | Units   | Description         | Note
------------|---------|---------------------|----------------------
$\delta$    | [rad]   | machine power angle |
$\omega$    | [p.u]   | machine speed       | Optionally read by a governor or a stabilizer component

#### Algebraic

Symbol  | Units  | Description                         | Note
--------|--------|-------------------------------------|-------------
$T_{e}$ | [p.u.] | electrical torque                   |
$I_r$   | [p.u.] | machine real injection current      | read by bus
$I_i$   | [p.u.] | machine imaginary injection current | read by bus

Note: All three can be expressed as function called by model equations. We add
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
\dot{\delta} &= (\omega - 1) \cdot \omega_0 \\
\dot{\omega} &= \frac{1}{2H}\left( \frac{P_{m} - D(\omega - 1)}{\omega}   - T_{e}\right)
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
    0 &= T_{e} - \frac{1}{\omega}\left( g E_p^2 - E_p \left[(gV_r - bV_i)\cos\delta + (bV_r + gV_i)\sin\delta \right]\right)\\
    0 &= I_r + gV_r - bV_i - E_p(g \cos\delta - b \sin\delta) \\
    0 &= I_i + bV_r + gV_i - E_p(b \cos\delta + g \sin\delta)
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

First compute injection currents from initial power injection power and bus
voltages:
```math
\begin{aligned}
I_r &= \frac{PV_r + QV_i}{V_r^2 + V_i^2} \\
I_i &= \frac{PV_i - QV_r}{V_r^2 + V_i^2}
\end{aligned} 
```

Next compute field winding voltage and machine angle:
```math
\begin{aligned}
E_r &= \frac{ g(I_r + gV_r - bV_i) + b (I_i + bV_r + gV_i)   }{g^2 + b^2} \\
E_i &= \frac{ -b(I_r + gV_r - bV_i) + g (I_i + bV_r + gV_i)   }{g^2 + b^2} \\
E_p &= \sqrt{E_r^2 + E_i^2} \\
\delta &= \arctan \dfrac{E_i}{E_r}
\end{aligned}
```

Set machine speed to the synchronous speed:
```math
\omega = 1
```

Now, we can compute electrical torque and set mechanical torque to be equal
to the electrical.
```math
\begin{aligned}
T_{elec} &= gE_p^2 - E_p \left[ (gV_r -  bV_i ) \cos\delta + (bV_r + gV_i )\sin\delta \right] \\
P_{mech} &= T_{elec} 
\end{aligned} 
```

With this, we initialize the machine at a steady state.
