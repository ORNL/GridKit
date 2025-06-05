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
  g   &=\dfrac{R}{R^2+(X)^2} \\
  b   &= -\dfrac{X}{R^2+(X)^2}\\
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
