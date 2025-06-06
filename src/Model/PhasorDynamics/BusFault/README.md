# Bus Fault

Represents an impedance fault at a bus. This device can exist in two states, on or off, controlled by the user. Following a state change, generally the solver needs to be reset as this is a discrete event.

## Model Parameters

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$R$  | [p.u.] | Fault resistance  | 
$X$  | [p.u.] | Fault reactance  | 
$U$ | [unitless] | Binary status $$\in \{0, 1\}$$ | Set by user to put fault on or off.

### Model Derived Parameters
``` math
\begin{aligned}
  G   &=\dfrac{R}{R^2+ X^2} \\
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
0 &= -I_{r} + U (-G V_{r} + B V_{i}) \\
0 &= -I_{i} + U (-B V_{r} - G V_{i})
\end{aligned}
```
