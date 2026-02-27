# Transmission Line Branch Model

Transmission lines and different types of transformers (traditional, Load
Tap-Changing transformers (LTC) and Phase Angle Regulators (PARs)) can be
modeled with a common branch model.

The most common circuit that is used to represent the transmission line model
is $`\pi`$ circuit as shown in Figure 1. The positive flow direction is into
buses. Commonly used convention is to define positive direction to be from
sending to receiving bus. We decide to use this symmetric convention because it
provides more flexibility for modeling.

<div align="center">
   <img align="center" src="../../../../docs/Figures/branch_phasor_dynamics.png">

  Figure 1: Transmission line $`\pi`$ equivalent circuit
</div>

## Model Parameters

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$R$  | [p.u.] | Branch series resistance  | 
$X$  | [p.u.] | Branch series reactance  | 
$G$  | [p.u.] | Branch shunt conductance  | 
$B$  | [p.u.] | Branch shunt susceptance  | 

### Model Derived Parameters
Note the difference between little-g and big-G, little-b, big-B in these equations.
``` math
\begin{aligned}
  g   &=\dfrac{R}{R^2 + X^2} \\
  b   &= -\dfrac{X}{R^2 + X^2}\\
\end{aligned}
```


## Model Variables

### Internal Variables

#### Differential
None.

#### Algebraic

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$I_{r1}$  | [p.u.] | Terminal current, real component, bus 1  | Read by bus
$I_{i1}$  | [p.u.] | Terminal current, imaginary component, bus 1  |  Read by bus
$I_{r2}$  | [p.u.] | Terminal current, real component, bus 2  | Read by bus
$I_{i2}$  | [p.u.] | Terminal current, imaginary component, bus 2  |  Read by bus


### External Variables

#### Differential
None.

#### Algebraic
Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$V_{r1}$  | [p.u.] | Terminal voltage, real component, bus 1 | owned by bus object
$V_{i1}$  | [p.u.] | Terminal voltage, imaginary component, bus 1 | owned by bus object
$V_{r2}$  | [p.u.] | Terminal voltage, real component, bus 2 | owned by bus object
$V_{i2}$  | [p.u.] | Terminal voltage, imaginary component, bus 2 | owned by bus object


## Model Equations

### Differential Equations
None.

### Algebraic Equations
``` math
\begin{aligned}
      0 &= - I_{r1} -\left(g + \dfrac{G}{2}\right) V_{r1} + \left(b + \dfrac{B}{2}\right) V_{i1} + g V_{r2} - b V_{i2}\\
      0 &= I_{i1} - \left(b + \dfrac{B}{2}\right) V_{r1} - \left(g + \dfrac{G}{2}\right) V_{i1} + b V_{r2} + g V_{i2}\\
      0 &= I_{r2} + g V_{r1} - b V_{i1} - \left(g + \dfrac{G}{2}\right) V_{r2} + \left(b + \dfrac{B}{2}\right) V_{i2}\\
      0 &= I_{i2} + b V_{r1} + g V_{i1} - \left(b + \dfrac{B}{2}\right) V_{r2} - \left(g + \dfrac{G}{2}\right) V_{i2}
\end{aligned}
```

# Model Outputs

2 outputs are model variables of the branch model: $I_r$ and $I_i$.

There are 3 calculated outputs for the branch model: current magnitude $|I|$, active power $P$ and reactive power $Q$. They are calculated as follows:
``` math
\begin{aligned}
      |I| &= \sqrt{(I_{r})^2 + (I_{i})^2} \\
      P &= V_{r} I_{r} + V_{i} I_{i}\\
      Q &= V_{i} I_{r} - V_{r} I_{i}
\end{aligned}
```
$|I|$ is the absolute value of the current phasor, independent of its direction or phase angle. The calculation above refelcts this as the geometric length.

In an AC electrical system, the complex power $S$ can be calculated from the voltage phasor $V$ and current phasor, $I$. By defenition, $S=VI^*$, where $I^*$ is the complex conjugate of $I$. Thus; 

``` math
\begin{aligned}
      S &= (V_r + V_i) (I_r - I_i)\\
      S &= V_r I_r - V_r I_i + V_i I_r + V_i I_i\\
      S &= (V_r I_r + V_i I_i) + (V_i I_r - V_r I_i)\\
      \Re({S}) &= V_r I_r + V_i I_i\\
      \Im({S}) &= V_i I_r - V_r I_i
\end{aligned}
```

$P$ is the real component of $S$. It is the power that actually does useful work. 

$Q$ is the imaginary component of $S$. It is the power that oscillates back and forth between capacitors and inductors. 

Positive $P_1$ values indicate that Bus 1 supplies active power into the branch. Negative values imply that Bus 1 is absorbing the active power. 
$\newline$
The same sign convention applies to reactive power: positive $Q_1$ corresponds to injection into the branch by Bus 1, and negative $Q_1$ corresponds to absorption from the branch.

The same output variables are computed for Bus 2 following the identical procedure used for Bus 1.

# Transformer Branch Model

**Note: Transformer model not yet implemented**

The branch model can be created by adding the ideal transformer in series with
the $`\pi`$ circuit as shown in Figure 2 where $`\tau`$ is a tap ratio
magnitude and $`\theta`$ is the phase shift angle and
$`N = \tau e^{j \theta}`$.

<div align="center">
   <img align="center" src="../../../../docs/Figures/transformer-branch.png">
   
   
  Figure 2: Branch equivalent circuit
</div>


The branch admitance matrix is then:

```math
\mathbf{Y}_{BR}=
\begin{bmatrix}
 -\left(g + jb + \dfrac{G+jB}{2} \right) \dfrac{1}{\tau^2} & (g + jb)\dfrac{1}{\tau e^{-j\theta}}\\
                                                           &                                     \\
     (g + jb)\dfrac{1}{\tau e^{j\theta}}                   & -\left(g + jb + \dfrac{G+jB}{2}\right)
\end{bmatrix}
```

### Branch contribution to residuals at adjacent busses

The currents entering adjacent buses are obtained in a similar manner as for
the $`\pi`$-model.
