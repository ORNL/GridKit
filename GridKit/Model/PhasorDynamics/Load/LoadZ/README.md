# LoadZ

Static constant-impedance load model. `LoadZ` has no solver-owned variables; it
computes terminal current contributions from the connected bus voltage and adds
them directly to the bus current-balance residuals.

## Model Parameters

Symbol | Units  | JSON | Description
-------|--------|------|------------
$R$    | [p.u.] | `R`  | Load resistance
$X$    | [p.u.] | `X`  | Load reactance

### Parameter Validation

None.

### Model Derived Parameters

```math
\begin{aligned}
G &= \frac{R}{R^2 + X^2} \\
B &= -\frac{X}{R^2 + X^2}
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

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|------
$V_{\mathrm{r}}$ | [p.u.] | Terminal voltage, real component      | Owned by connected bus
$V_{\mathrm{i}}$ | [p.u.] | Terminal voltage, imaginary component | Owned by connected bus

## Wiring

Port  | Type | Description
------|------|------------
`bus` | Bus  | Connected bus that owns terminal voltage variables and current-balance residuals

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Bus Current-Balance Contributions

Let $I_{\mathrm{r}}^{\mathrm{LoadZ}}$ and
$I_{\mathrm{i}}^{\mathrm{LoadZ}}$ denote the model contributions to the real
and imaginary current-balance residuals of the connected bus. Positive current
is oriented entering the bus.

```math
\begin{aligned}
I_{\mathrm{r}}^{\mathrm{LoadZ}}
  &= -(G V_{\mathrm{r}} - B V_{\mathrm{i}}) \\
I_{\mathrm{i}}^{\mathrm{LoadZ}}
  &= -(B V_{\mathrm{r}} + G V_{\mathrm{i}})
\end{aligned}
```

These contributions are accumulated directly into the bus-owned residuals.

## Initialization

None.

## Monitors

Name | Units  | Description                                  | Note
-----|--------|----------------------------------------------|------
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
