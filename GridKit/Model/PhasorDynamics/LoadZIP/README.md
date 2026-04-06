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
0 &= I_{r} + (P_{0} V_{r} + Q_{0} V_{i}) \left(\frac{1}{V_0^2} (1 - \alpha_I - \alpha_P) + \frac{1}{V_0 \sqrt{V_r^2+V_i^2}} \alpha_I + \frac{1}{V_r^2+V_I^2} \alpha_P\right) \\
0 &= I_{i} + (P_{0} V_{i} - Q_{0} V_{r}) \left(\frac{1}{V_0^2} (1 - \alpha_I - \alpha_P) + \frac{1}{V_0 \sqrt{V_r^2+V_i^2}} \alpha_I + \frac{1}{V_r^2+V_I^2} \alpha_P\right)
\end{aligned}
```

### Initialization Procedure
Use the algebraic equations to solve for $I_{r}$ and $I_{i}$.

### Model Outputs
Real and imaginary values of the load current are the variables 
$I_{r}$ and $I_{i}$.

Current is oriented leaving the load (i.e. entering the bus).

Current magnitude $I_{m}$ is the phasor magnitude of the current.
``` math
\begin{aligned}
      I_{m} &= \sqrt{(I_{r})^2 + (I_{i})^2}
\end{aligned}
```

Active and reactive power ($P$ and $Q$)
are the real and imaginary parts of the complex power,
where the complex power is defined as $S=VI^{\ast}=(V_r + j V_i)(I_r - jI_i)$
``` math
\begin{aligned}
      P &= V_{r} I_{r} + V_{i} I_{i}\\
      Q &= V_{i} I_{r} - V_{r} I_{i}
\end{aligned}
```


### Note on Derivation
The origin of the algebraic equations is easier to understand in complex form:
``` math
\begin{aligned}
S &= S_z  \left(\frac{|V|}{V_0}\right)^2 + S_I \left(\frac{|V|}{V_0}\right) + S_P \\
S &= V I^*
\end{aligned}
```