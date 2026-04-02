# LoadZIP Model

This represents a static load model with impedance (Z), current (I), and
power (P) components.

## Model Parameters

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$P_0$  | [p.u.] | Load nominal real power  | 
$Q_0$  | [p.u.] | Load nominal reactive power  | 
$V_0$  | [p.u.] | Load nominal voltage | 
$\alpha_I$  | [unitless] | Fraction of load to be represented as constant current | 
$\alpha_P$  | [unitless] | Fraction of load to be represented as constant power | 



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
0 &= I_{r} - (P_{0} V_{r} + Q_{0} V_{i}) \left(\frac{1}{V_0^2} (1 - \alpha_I - \alpha_P) + \frac{1}{V_0 \sqrt{V_r^2+V_i^2}} \alpha_I + \frac{1}{V_r^2+V_I^2} \alpha_P\right) \\
0 &= I_{i} - (P_{0} V_{i} - Q_{0} V_{r}) \left(\frac{1}{V_0^2} (1 - \alpha_I - \alpha_P) + \frac{1}{V_0 \sqrt{V_r^2+V_i^2}} \alpha_I + \frac{1}{V_r^2+V_I^2} \alpha_P\right)
\end{aligned}
```


### Note on Derivation
The origin of the algebraic equations is easier to understand in complex form:
``` math
\begin{aligned}
S &= S_z  \left(\frac{|V|}{V_0}\right)^2 + S_I \left(\frac{|V|}{V_0}\\right) + S_P
S &= V I^*
\end{aligned}
```