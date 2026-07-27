# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` has no solver-owned variables; it computes terminal
current contributions from the connected bus voltage and adds them directly to
the bus current-balance residuals.

## Model Parameters

Symbol          | Units  | JSON     | Description                    | Typical Value
----------------|--------|----------|--------------------------------|--------------
$P_\text{nom}$  | [p.u.] | `Pnom`   | Nominal consumed real power    | 0.0
$Q_\text{nom}$  | [p.u.] | `Qnom`   | Nominal consumed reactive power | 0.0
$\alpha_I$      | [-]    | `alphaI` | Constant current load fraction | 0.0
$\alpha_P$      | [-]    | `alphaP` | Constant power load fraction   | 0.0

### Parameter Validation

None.

### Model Derived Parameters

$V_\text{nom}$ is the initial voltage magnitude of the respective bus.

```math
\begin{aligned}
G &= \frac{P_{\mathrm{nom}}}{V_{\mathrm{nom}}^2} \\
B &= \frac{Q_{\mathrm{nom}}}{V_{\mathrm{nom}}^2} \\
\alpha_{Z} &= 1 - \alpha_{I} - \alpha_{P}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$V_{\mathrm{r}}$ | [p.u.] | Terminal voltage, real component      | Owned by connected bus
$V_{\mathrm{i}}$ | [p.u.] | Terminal voltage, imaginary component | Owned by connected bus

## Wiring

Port  | Type | Description
------|------|------------
`bus` | Bus  | Connected bus that owns terminal voltage variables and current-balance residuals

## Model Equations

Let $V = \sqrt{V_{\mathrm{r}}^2 + V_{\mathrm{i}}^2}$.

### Differential Equations

None.

### Algebraic Equations

None.

### Bus Current-Balance Contributions

Let $I_{\mathrm{r}}^{\mathrm{LoadZIP}}$ and
$I_{\mathrm{i}}^{\mathrm{LoadZIP}}$ denote the model contributions to the real
and imaginary current-balance residuals of the connected bus. Positive current
is oriented entering the bus.

```math
\begin{aligned}
I_{\mathrm{r}}^{\mathrm{LoadZIP}}
  &= -(G V_{\mathrm{r}} + B V_{\mathrm{i}})
\left[
\alpha_{Z}
+ \alpha_{I} \frac{V_{\mathrm{nom}}}{V}
+ \alpha_{P} \frac{V_{\mathrm{nom}}^2}{V^2}
\right] \\
I_{\mathrm{i}}^{\mathrm{LoadZIP}}
  &= -(G V_{\mathrm{i}} - B V_{\mathrm{r}})
\left[
\alpha_{Z}
+ \alpha_{I} \frac{V_{\mathrm{nom}}}{V}
+ \alpha_{P} \frac{V_{\mathrm{nom}}^2}{V^2}
\right]
\end{aligned}
```

These contributions are accumulated directly into the bus-owned residuals.

## Initialization

```math
V_\text{nom} \leftarrow \sqrt{V_\mathrm{r}^2 + V_\mathrm{i}^2}
```

The nominal-voltage anchor and derived admittance parameters are recomputed
from the initialized bus voltage. The model has no internal state to initialize.

## Monitors

Name | Units  | Description                                  | Note
-----|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
