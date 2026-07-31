# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` owns terminal current states and adds their current
contribution to the connected bus residual.

## Model Parameters

Symbol          | Units  | JSON     | Description                    | Typical Value
----------------|--------|----------|--------------------------------|--------------
$P_\text{nom}$  | [p.u.] | `Pnom`   | Legacy nominal real-power fallback | 0.0
$Q_\text{nom}$  | [p.u.] | `Qnom`   | Legacy nominal reactive-power fallback | 0.0
$V_\text{nom}$  | [p.u.] | `Vnom`   | Nominal voltage magnitude      | 1.0
$\alpha_I$      | [-]    | `alphaI` | Constant current load fraction | 0.0
$\alpha_P$      | [-]    | `alphaP` | Constant power load fraction   | 0.0

### Parameter Validation

None.

### Model Derived Parameters

When the legacy nominal-power fallback is used,

```math
\begin{aligned}
G &= \dfrac{P_\mathrm{nom}}{V_\mathrm{nom}^2} \\
B &= -\dfrac{Q_\mathrm{nom}}{V_\mathrm{nom}^2} \\
\alpha_Z &= 1 - \alpha_I - \alpha_P
\end{aligned}
```

## Model Inputs

Input    | Units  | Description                                      | Default
---------|--------|--------------------------------------------------|--------
`p`      | [p.u.] | Initial terminal active-power injection          | Derived from `Pnom`
`q`      | [p.u.] | Initial terminal reactive-power injection        | Derived from `Qnom`
`online` | [-]    | Connection status; zero is disconnected          | 1

Positive `p` and `q` inject power into the bus. They are sampled during
initialization. Any nonzero `online` value connects the load; status is read
when the network contribution and monitors are evaluated. Disconnecting the
load does not remove its internal algebraic equations. Change inputs only
while the solve is stopped.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$I_r$  | [p.u.] | Terminal current, real component         | Added to connected bus residual
$I_i$  | [p.u.] | Terminal current, imaginary component    | Added to connected bus residual

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$V_r$  | [p.u.] | Terminal voltage, real component         | Owned by connected bus
$V_i$  | [p.u.] | Terminal voltage, imaginary component    | Owned by connected bus

## Model Equations

Let $V = \sqrt{V_r^2 + V_i^2}$.

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= I_r + (G V_r - B V_i)
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right] \\
0 &= I_i + (G V_i + B V_r)
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right]
\end{aligned}
```

### External Equations

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i
\end{aligned}
```

## Initialization

Initialization solves the algebraic current states from the connected bus
voltage. Let $V_{m0} = \sqrt{V_{r0}^2 + V_{i0}^2}$.

When terminal inputs $p$ and $q$ are attached, the model first chooses $G$ and
$B$ so those injections are reproduced at the initialized voltage:

```math
\begin{aligned}
z_0 &= \alpha_Z + \alpha_I \frac{V_\text{nom}}{V_{m0}}
     + \alpha_P \frac{V_\text{nom}^2}{V_{m0}^2} \\
G &= -\frac{p}{V_{m0}^2 z_0}, \qquad
B = -\frac{q}{V_{m0}^2 z_0}.
\end{aligned}
```

```math
\begin{aligned}
    I_r &\leftarrow -G V_r + B V_i \\
    I_i &\leftarrow -G V_i - B V_r
\end{aligned}
```

The derivative vector entries initialize to zero.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
