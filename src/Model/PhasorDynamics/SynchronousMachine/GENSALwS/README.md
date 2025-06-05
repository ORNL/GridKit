# GENSAL 

This synchronous machine model is 5th order and is specifically designed for salient pole machines. It is a standard model used in phasor-domain industry stability studies.
See the [General Synchronous Machine Model](../README.md) for general synchronous machine information.

Notes:
- $X''_{q}=X''_{d}$  (no subtransient saliency)
- $X''_{d}$ does not saturate
- Only d-axis affected by saturation
- $X_{q}=X'_{q}$
- $T'_{q0}$ is neglected

## Block Diagram
<div align="center">
   <img align="center" src="../../../../../docs/Figures/GENSAL.JPG">
   
   
  Figure 2: GENSAL. Figure courtesy of 
  [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol      | Units   | Description                     | Typical Value | Note
------------|---------|---------------------------------|---------------| ------
$\omega_0$  | [rad/s] | synchronous frequency           | $2\pi60$
$H$         | [s]     | rotor inertia                   | 3
$D$         | [p.u.]  | damping coefficient             | 0
$R_a$       | [p.u.]  | winding resistance              | 0
$X_{\ell}$   | [p.u.] | Stator leakage reactance | 0.15 | 
$X_{d}$   | [p.u.] | Direct axis synchronous reactance | 2.1 | 
$X'_{d}$   | [p.u.] | Direct axis transient reactance | 0.2 | 
$X''_{d}$   | [p.u.] | Direct axis sub-transient reactance | 0.18 | 
$X_{q}$   | [p.u.] | Quadrature axis synchronous reactance | 0.5 | 
$T'_{d0}$   | [s] | Open circuit direct axis transient time const. | 7 | 
$T''_{d0}$   | [s] | Open circuit direct axis sub-transient time const. | 0.04 | 
$T''_{q0}$   | [s] | Open circuit quadrature axis sub-transient time const. | 0.05 | 
$S_{10}$   | [p.u.] | Saturation factor at 1.0 pu flux | 0 | 
$S_{12}$   | [p.u.] | Saturation factor at 1.2 pu flux | 0 | 

### Model Derived Parameters
``` math
\begin{aligned}
  G   &=\dfrac{R_a}{R_a^2+(X''_q)^2}&
  B   &= -\dfrac{X''_q}{R_a^2+(X''_q)^2}\\
  S_A &= \dfrac{1.2\sqrt{S_{10}/S_{12}} +1}{\sqrt{S_{10}/S_{12}} +1} & 
  S_B &= \dfrac{1.2\sqrt{S_{10}/S_{12}} -1}{\sqrt{S_{10}/S_{12}} -1} \\
  X_{d1} &= X_d-X_d'       \\
  X_{d2} &= X_d'-X_\ell & X_{q2} = X_q-X''_q \\
  X_{d3} &= (X_d'-X_d'')/X_{d2}^2  \\
  X_{d5} &= (X_d''-X_\ell)/X_{d2}  \\
  X_{qd} &= (X_q-X_\ell)/(X_d-X_\ell)
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$\delta$  | [rad] | Machine internal rotor angle |
$\omega$  | [p.u.] | Machine relative speed | Optionally read by governor or stabilizer component
$\psi'_d$ | [p.u.] | Direct axis subtransient flux | 
$\psi''_q$ | [p.u.] | Quadrature axis subtransient flux | 
$E'_q$ | [p.u.] | Quadrature axis subtransient flux | 

#### Algebraic
Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$V_d$  | [p.u.] | Machine internal voltage, d-axis | 
$V_q$  | [p.u.] | Machine internal voltage, q-axis | 
$I_d$  | [p.u.] | Terminal current, d-axis  | 
$I_q$  | [p.u.] | Terminal current, q-axis   | 
$I_r$  | [p.u.] | Terminal current, real component on network reference frame  | Read by bus and optionally by controllers
$I_i$  | [p.u.] | Terminal current, imaginary component on network reference frame  |  Read by bus and optionally by controllers
$\psi''_d$  | [p.u.] | Total d-axis subtransient flux
$T_{elec}$   | [p.u.] |  Electrical torque

### External Variables

#### Differential
None.

#### Algebraic
Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$V_r$  | [p.u.] | Terminal voltage, real component on network reference frame | owned by bus object
$V_i$  | [p.u.] | Terminal voltage, imaginary component on network reference frame | owned by bus object
$P_{mech}$   | [p.u.] | Mechanical power from the prime mover | Owned by governor, constant if no governor is connected to the machine
$E_{fd}$    | [p.u.] | Field winding voltage from the excitation system  | Owned by exciter, constant if no exciter is connected to the machine

## Model Equations

### Differential Equations
``` math
\begin{aligned}
  \dot\delta      &= \omega\cdot\omega_0 \\
  \dot\omega      &= \dfrac{1}{2H}\left(\dfrac{P_{mech}-D\omega}{1+\omega}
                   - T_{elec}\right)\\
  \dot{\psi}'_{d} &= \dfrac{1}{T''_{d0}}(E'_{q}-\psi'_{d}-X_{d2}I_{d})\\
  \dot{\psi}''_{q} &= \dfrac{1}{T''_{q0}}(-\psi''_{q}-X_{q2}I_{q})\\
  \dot{E}'_{q} &= \dfrac{1}{T'_{d0}}
    \left(
      E_{fd}-E'_{q}-X_{d1}
      (I_{d}+X_{d3}(E'_{q}-\psi'_{d}-X_{d2}I_{d}))
      - S_B (E'_q-S_A)^2
    \right)\\
\end{aligned}
```

### Algebraic Equations
Note that for implementation purposes, some of these equations may be simplified into functions and the internal variables eliminated. Nevertheless, for modeling clarity and conformance to typical practice, the full equations are given here.
``` math
\begin{aligned}
  0 &= -V_{d} -\psi''_{q}(1+\omega)\\
  0 &= -V_{q}  +\psi''_{d}(1+\omega)\\
 0 &= -I_d + I_r \sin(\delta) - I_i \cos(\delta) \\
 0 &= -I_q + I_r \cos(\delta) + I_i \sin(\delta) \\
 0 &= -I_r + G (V_d \sin(\delta) + V_q \cos(\delta) - V_r) - B (V_d \cos(\delta) + V_q \sin(\delta) - V_i) \\
 0 &= -I_i + B (V_d \sin(\delta) + V_q \cos(\delta) - V_r) + G (V_d \cos(\delta) + V_q \sin(\delta) - V_i) \\
  0 &= -\psi''_{d} +E'_{q}X_{d5} + \psi'_{d}X_{d4}\\
  0 &= -T_{elec} +(\psi''_{d} - I_dX_d'')I_q-(\psi''_{q} - I_qX_d'')I_d
\end{aligned}
```
