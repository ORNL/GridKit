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
  S_A    &= \min\left(\dfrac{1.2\sqrt{S_{10}/S_{12}} +1}{\sqrt{S_{10}/S_{12}} +1},
                       \dfrac{1.2\sqrt{S_{10}/S_{12}} -1}{\sqrt{S_{10}/S_{12}} -1}\right) &
  S_B    &= \dfrac{S_{12}}{(S_A-1.2)^2} \\
  X_{d1} &= X_d-X_d'                 & X_{q2} &= X_q-X_d'' \\
  X_{d2} &= X_d'-X_\ell              & X_{d3} &= (X_d'-X_d'')/X_{d2}^2 \\
  X_{d4} &= (X_d'-X_d'')/X_{d2}      & X_{d5} &= (X_d''-X_\ell)/X_{d2} \\
  f_\mathrm{base} &= f_\mathrm{sys} &
  S_\mathrm{mach,VA} &= 10^6 S_\mathrm{mach}
\end{aligned}
```

When $S_{12}=0$, $S_A=S_B=0$.

System bases are taken from the system at initialization.

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------------
`bus`   | Bus    | Known   | Terminal bus voltage and current-balance residuals
`pmech` | Input  | Unknown | System-base mechanical-power input; converted to machine base internally and held constant when unconnected
`efd`   | Input  | Unknown | Machine-base field-voltage input; held constant when unconnected
`speed` | Output | Known   | Machine speed-deviation output

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
$k_{sat}$  | [p.u.] | Saturation factor                 |
$V_d$      | [p.u.] | Machine internal voltage, d-axis  |
$V_q$      | [p.u.] | Machine internal voltage, q-axis  |
$T_e$      | [p.u.] | Electrical torque                 |
$I_d$      | [p.u.] | Terminal current, d-axis          |
$I_q$      | [p.u.] | Terminal current, q-axis          |
$I_r$      | [p.u.] | Terminal current, real component on network reference frame | Machine base; converted to system base for the bus and monitors
$I_i$      | [p.u.] | Terminal current, imaginary component on network reference frame | Machine base; converted to system base for the bus and monitors

### External Variables

#### Differential
None.

#### Algebraic
Symbol   | Units  | Description                                             | Note
---------|--------|---------------------------------------------------------| ------
$V_r$    | [p.u.] | Terminal voltage, real component on network reference frame      | owned by bus object
$V_i$    | [p.u.] | Terminal voltage, imaginary component on network reference frame | owned by bus object
$P_m$    | [p.u.] | Mechanical power from the prime mover                   | System-base signal; converted to machine base internally and held constant if unconnected
$E_{fd}$ | [p.u.] | Field winding voltage from the excitation system        | Machine-base signal; held constant if unconnected

## Model Equations

### Internal Equations

#### Differential

``` math
\begin{aligned}
  \dot\delta       &= \omega \cdot 2\pi f_\mathrm{base} \\
  \dot\omega       &= \dfrac{1}{2H}\left(\dfrac{P_m-D\omega}{1+\omega}
                    - T_e\right)\\
  \dot{E}'_q       &= \dfrac{1}{T'_{d0}}
    \left(
      E_{fd}-E'_q-X_{d1}
      (I_d+X_{d3}(E'_q-\psi'_d-X_{d2}I_d))
      -E'_q k_{sat}
    \right)\\
  \dot{\psi}'_d    &= \dfrac{1}{T''_{d0}}(E'_q-\psi'_d-X_{d2}I_d)\\
  \dot{\psi}''_q   &= \dfrac{1}{T''_{q0}}(-\psi''_q-X_{q2}I_q)
\end{aligned}
```

#### Algebraic

``` math
\begin{aligned}
  0 &= -\psi''_d + E'_qX_{d5}+\psi'_dX_{d4}\\
  0 &= -k_{sat} + S_B q(E'_q-S_A)\\
  0 &= -V_d -\psi''_q(1+\omega)\\
  0 &= -V_q +\psi''_d(1+\omega)\\
  0 &= -T_e +(\psi''_d-I_dX_d'')I_q-(\psi''_q-I_qX_d'')I_d\\
  0 &= -I_d + I_r \sin(\delta) - I_i \cos(\delta) \\
  0 &= -I_q + I_r \cos(\delta) + I_i \sin(\delta) \\
  0 &= -I_r + G (V_d \sin(\delta) + V_q \cos(\delta) - V_r) - B (-V_d \cos(\delta) + V_q \sin(\delta) - V_i) \\
  0 &= -I_i + B (V_d \sin(\delta) + V_q \cos(\delta) - V_r) + G (-V_d \cos(\delta) + V_q \sin(\delta) - V_i)
\end{aligned}
```

CommonMath defines the primitive
[quadratic ramp](../../../../CommonMath.md#quadratic-ramp) $q$.

### External Equations

The machine-base terminal currents are converted to system base and added to
the connected bus residuals. Here $I_r^{\mathrm{mach}}\equiv I_r$ and
$I_i^{\mathrm{mach}}\equiv I_i$ denote the internal machine-base currents, and
$S_\mathrm{sys,VA}$ is the system power base in volt-amperes:

```math
\begin{aligned}
I_r^{\mathrm{bus}}
  &\leftarrow I_r^{\mathrm{bus}}
  + \dfrac{S_\mathrm{mach,VA}}{S_\mathrm{sys,VA}} I_r^{\mathrm{mach}} \\
I_i^{\mathrm{bus}}
  &\leftarrow I_i^{\mathrm{bus}}
  + \dfrac{S_\mathrm{mach,VA}}{S_\mathrm{sys,VA}} I_i^{\mathrm{mach}}.
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
  k_{sat} &= S_B q(E'_q-S_A)\\
  T_e &= (\psi''_d-I_dX_d'')I_q-(\psi''_q-I_qX_d'')I_d\\
  P_m &= T_e\\
  E_{fd} &= E'_q+X_{d1}(I_d+X_{d3}(E'_q-\psi'_d-X_{d2}I_d))+E'_q k_{sat}
\end{aligned}
```

## Monitors

Monitor | Units | Description                                                        | Note
--------|-------|--------------------------------------------------------------------|------
`ir`     | [p.u.] | Terminal current, real component $I_r$ in the network frame         | Oriented leaving the machine; system base
`ii`     | [p.u.] | Terminal current, imaginary component $I_i$ in the network frame    | Oriented leaving the machine; system base
`p`      | [p.u.] | Active power $P=V_rI_r+V_iI_i$                                     | Oriented leaving the machine; system base
`q`      | [p.u.] | Reactive power $Q=V_iI_r-V_rI_i$                                   | Oriented leaving the machine; system base
`delta`  | [rad]  | Machine internal rotor angle $\delta$                               |
`omega`  | [p.u.] | Machine speed deviation $\omega$                                   | $\omega=0$ at synchronous speed
`speed`  | [p.u.] | Per-unit machine speed                                              | $1+\omega$
`Eqp`    | [p.u.] | Quadrature-axis transient flux $E'_q$                               | Machine base
`psidp`  | [p.u.] | Direct-axis transient flux $\psi'_d$                               | Machine base
`psiqpp` | [p.u.] | Total q-axis subtransient flux $\psi''_q$                           | Machine base
`psidpp` | [p.u.] | Total d-axis subtransient flux $\psi''_d$                           | Machine base
`vd`     | [p.u.] | Machine internal voltage, d-axis $V_d$                              | Machine base
`vq`     | [p.u.] | Machine internal voltage, q-axis $V_q$                              | Machine base
`te`     | [p.u.] | Electrical torque $T_e$                                             | Machine base
`id`     | [p.u.] | Terminal current, d-axis $I_d$                                      | Machine base
`iq`     | [p.u.] | Terminal current, q-axis $I_q$                                      | Machine base
