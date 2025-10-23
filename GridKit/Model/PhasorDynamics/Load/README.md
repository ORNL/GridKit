# Load Model

Load modeling is one of the more complex aspects of power system dynamics.
The simplest model, which is used for this challenge problem, is to model
the load as a complex shunt impedance $$ R + jX $$.


## Model Parameters

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$R$  | [p.u.] | Load resistance  | 
$X$  | [p.u.] | Load reactance  | 


### Model Derived Parameters
``` math
\begin{aligned}
  G   &=\dfrac{R}{R^2 + X^2} \\
  B   &= -\dfrac{X}{R^2 + X^2}\\
\end{aligned}
```


## Model Variables

### Internal Variables

#### Differential
None.

#### Algebraic

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$I_r$  | [p.u.] | Terminal current, real component  | Read by bus
$I_i$  | [p.u.] | Terminal current, imaginary component  |  Read by bus


### External Variables

#### Differential
None.

#### Algebraic
Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$V_r$  | [p.u.] | Terminal voltage, real component | owned by bus object
$V_i$  | [p.u.] | Terminal voltage, imaginary component | owned by bus object


## Model Equations

### Differential Equations
None.

### Algebraic Equations
``` math
\begin{aligned}
0 &= I_{r} +G V_{r} - B V_{i} \\
0 &= I_{i} +B V_{r} + G V_{i}
\end{aligned}
```
