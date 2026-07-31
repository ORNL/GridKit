# GENSAL

This synchronous machine model is 5th order and is specifically designed for
salient-pole machines. It is a standard model used in phasor-domain industry
stability studies.
See the [General Synchronous Machine Model](../README.md) for general
synchronous machine information.

Notes:
- $X_q''=X_d''$ (no subtransient saliency)
- $X_q=X_q'$
- $T'_{q0}$ is neglected
- Only d-axis affected by saturation

## Block Diagram
![](../../../../../docs/Figures/GENSAL.JPG)

Figure 2: GENSAL. Figure courtesy of
[PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol      | Units   | Description                                | Typical Value | Note
------------|---------|--------------------------------------------|---------------| ------
$P_0$       | [p.u.]  | Initial active power injection             | 1.0 |
$Q_0$       | [p.u.]  | Initial reactive power injection           | 0.0 |
$H$         | [s]     | rotor inertia                              | 3
$D$         | [p.u.]  | damping coefficient                        | 0
$R_a$       | [p.u.]  | winding resistance                         | 0
$T'_{d0}$   | [s]     | Open circuit direct axis transient time const. | 7 |
$T''_{d0}$  | [s]     | Open circuit direct axis sub-transient time const. | 0.04 |
$T''_{q0}$  | [s]     | Open circuit quadrature axis sub-transient time const. | 0.05 |
$X_d$       | [p.u.]  | Direct axis synchronous reactance          | 2.1 |
$X'_d$      | [p.u.]  | Direct axis transient reactance            | 0.2 |
$X''_d$     | [p.u.]  | Direct axis sub-transient reactance        | 0.18 |
$X_q$       | [p.u.]  | Quadrature axis synchronous reactance      | 0.5 |
$X_{\ell}$  | [p.u.]  | Stator leakage reactance                   | 0.15 |
$S_{10}$    | [p.u.]  | Saturation factor at 1.0 pu flux           | 0 |
$S_{12}$    | [p.u.]  | Saturation factor at 1.2 pu flux           | 0 |
$S_\mathrm{mach}$ | [MVA] | Machine power base                  | 100 |

### Model Derived Parameters
``` math
\begin{aligned}
  G      &=  \dfrac{R_a}{R_a^2+(X_d'')^2} &
  B      &= -\dfrac{X_d''}{R_a^2+(X_d'')^2}\\
  S_A    &= \dfrac{1.2\sqrt{S_{10}/S_{12}} +1}{\sqrt{S_{10}/S_{12}} +1} &
  S_B    &= \dfrac{1.2\sqrt{S_{10}/S_{12}} -1}{\sqrt{S_{10}/S_{12}} -1} \\
  X_{d1} &= X_d-X_d'                 & X_{q2} &= X_q-X_d'' \\
  X_{d2} &= X_d'-X_\ell              & X_{d3} &= (X_d'-X_d'')/X_{d2}^2 \\
  X_{d4} &= (X_d'-X_d'')/X_{d2}      & X_{d5} &= (X_d''-X_\ell)/X_{d2} \\
  f_\mathrm{base} &= f_\mathrm{sys} &
  S_\mathrm{mach,VA} &= 10^6 S_\mathrm{mach}
\end{aligned}
```

System bases are taken from the system at initialization.

## Model Variables

### Internal Variables

#### Differential

Symbol       | Units  | Description                       | Note
-------------|--------|-----------------------------------|-------
$\delta$     | [rad]  | Machine internal rotor angle      |
$\omega$     | [p.u.] | Machine speed deviation           | Optionally read by governor or stabilizer component
$E'_q$       | [p.u.] | Quadrature axis transient flux    |
$\psi'_d$    | [p.u.] | Direct axis transient flux        |
$\psi''_q$   | [p.u.] | Total q-axis subtransient flux    |

#### Algebraic
Symbol      | Units  | Description                       | Note
------------|--------|-----------------------------------| ------
$\psi''_d$ | [p.u.] | Total d-axis subtransient flux    |
$k_{sat}$  | [p.u.] | Additive saturation signal        |
$V_d$      | [p.u.] | Machine internal voltage, d-axis  |
$V_q$      | [p.u.] | Machine internal voltage, q-axis  |
$T_e$      | [p.u.] | Electrical torque                 |
$I_d$      | [p.u.] | Terminal current, d-axis          |
$I_q$      | [p.u.] | Terminal current, q-axis          |
$I_r$      | [p.u.] | Terminal current, real component on network reference frame | Read by bus and optionally by controllers
$I_i$      | [p.u.] | Terminal current, imaginary component on network reference frame | Read by bus and optionally by controllers

### External Variables

#### Differential
None.

#### Algebraic
Symbol   | Units  | Description                                             | Note
---------|--------|---------------------------------------------------------| ------
$V_r$    | [p.u.] | Terminal voltage, real component on network reference frame      | owned by bus object
$V_i$    | [p.u.] | Terminal voltage, imaginary component on network reference frame | owned by bus object
$P_m$    | [p.u.] | Mechanical power from the prime mover                   | Owned by governor, constant if no governor is connected to the machine
$E_{fd}$ | [p.u.] | Field winding voltage from the excitation system        | Owned by exciter, constant if no exciter is connected to the machine

## Model Equations

### Differential Equations
``` math
\begin{aligned}
  \dot\delta       &= \omega \cdot 2\pi f_\mathrm{base} \\
  \dot\omega       &= \dfrac{1}{2H}\left(\dfrac{P_m-D\omega}{1+\omega}
                    - T_e\right)\\
  \dot{E}'_q       &= \dfrac{1}{T'_{d0}}
    \left(
      E_{fd}-E'_q-X_{d1}
      (I_d+X_{d3}(E'_q-\psi'_d-X_{d2}I_d))
      -k_{sat}
    \right)\\
  \dot{\psi}'_d    &= \dfrac{1}{T''_{d0}}(E'_q-\psi'_d-X_{d2}I_d)\\
  \dot{\psi}''_q   &= \dfrac{1}{T''_{q0}}(-\psi''_q-X_{q2}I_q)
\end{aligned}
```

### Algebraic Equations
``` math
\begin{aligned}
  0 &= -\psi''_d + E'_qX_{d5}+\psi'_dX_{d4}\\
  0 &= -k_{sat} + S_B(E'_q-S_A)^2\sigma(E'_q-S_A)\\
  0 &= -V_d -\psi''_q(1+\omega)\\
  0 &= -V_q +\psi''_d(1+\omega)\\
  0 &= -T_e +(\psi''_d-I_dX_d'')I_q-(\psi''_q-I_qX_d'')I_d\\
  0 &= -I_d + I_r \sin(\delta) - I_i \cos(\delta) \\
  0 &= -I_q + I_r \cos(\delta) + I_i \sin(\delta) \\
  0 &= -I_r + G (V_d \sin(\delta) + V_q \cos(\delta) - V_r) - B (-V_d \cos(\delta) + V_q \sin(\delta) - V_i) \\
  0 &= -I_i + B (V_d \sin(\delta) + V_q \cos(\delta) - V_r) + G (-V_d \cos(\delta) + V_q \sin(\delta) - V_i)
\end{aligned}
```

## Initialization

Using the power-flow solution, initial currents are calculated from active and
reactive power injection. The remaining variables are initialized from the
steady-state GENSAL equations.

``` math
\begin{aligned}
  \omega &= 0 \\
  \delta &= \text{arg}\left[V_r+jV_i+(R_a+jX_q)(I_r+jI_i)\right]\\
  I_d &= I_r\sin(\delta)-I_i\cos(\delta)\\
  I_q &= I_r\cos(\delta)+I_i\sin(\delta)\\
  \psi''_q &= -X_{q2}I_q\\
  V_d &= -\psi''_q\\
  V_q &= V_r\cos(\delta)+V_i\sin(\delta)+X_d''I_d+R_aI_q\\
  \psi''_d &= V_q\\
  \psi'_d &= \psi''_d-(X_d''-X_\ell)I_d\\
  E'_q &= \psi'_d+X_{d2}I_d\\
  k_{sat} &= S_B(E'_q-S_A)^2\sigma(E'_q-S_A)\\
  T_e &= (\psi''_d-I_dX_d'')I_q-(\psi''_q-I_qX_d'')I_d\\
  P_m &= T_e\\
  E_{fd} &= E'_q+X_{d1}(I_d+X_{d3}(E'_q-\psi'_d-X_{d2}I_d))+k_{sat}
\end{aligned}
```

## Model Outputs

Symbol     | Units  | Description                       | Note
-----------|--------|-----------------------------------|------
$I_r$      | [p.u.] | Terminal current, real component on network reference frame | Oriented leaving the machine, system base
$I_i$      | [p.u.] | Terminal current, imaginary component on network reference frame | Oriented leaving the machine, system base
$P$        | [p.u.] | Active power, $V_rI_r+V_iI_i$     | Oriented leaving the machine, system base
$Q$        | [p.u.] | Reactive power, $V_iI_r-V_rI_i$   | Oriented leaving the machine, system base
$\delta$   | [rad]  | Machine internal rotor angle      |
$\omega$   | [p.u.] | Machine speed deviation           | $\omega=0$ at synchronous speed
$\text{speed}$ | [p.u.] | Per-unit machine speed            | $1+\omega$
$E'_q$     | [p.u.] | Quadrature axis transient flux    | Machine base
$\psi'_d$  | [p.u.] | Direct axis transient flux        | Machine base
$\psi''_q$ | [p.u.] | Total q-axis subtransient flux    | Machine base
$\psi''_d$ | [p.u.] | Total d-axis subtransient flux    | Machine base
$V_d$      | [p.u.] | Machine internal voltage, d-axis  | Machine base
$V_q$      | [p.u.] | Machine internal voltage, q-axis  | Machine base
$T_e$      | [p.u.] | Electrical torque                 | Machine base
$I_d$      | [p.u.] | Terminal current, d-axis          | Machine base
$I_q$      | [p.u.] | Terminal current, q-axis          | Machine base
