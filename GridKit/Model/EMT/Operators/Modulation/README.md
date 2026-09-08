# Modulation Model

`Modulation` normalizes three-phase converter voltage commands by the DC-link
voltage for sinusoidal [PWM](../../Component/Controller/PWM/README.md).

## Block Diagram

None.

## Model Parameters

None.

### Parameter Validation

None.

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `u` | Input | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^3$
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC-link voltage | $v_{\mathrm{dc}} > 0$
$\mathbf{m}$ | `m` | Output | [-] | Modulation command | $\mathbf{m} \in \mathbb{R}^3$

Vectors use $(a,b,c)$ order. All inputs must be connected and finite.
The upstream voltage controller limits the command to the available DC voltage.

## Submodels

None.

### Submodel Validation

None.

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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^3$
$v_{\mathrm{dc}}$ | [V] | DC-link voltage | Positive

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\mathbf{m} \leftarrow \dfrac{2}{v_{\mathrm{dc}}}\mathbf{u}
```

The output is an algebraic expression without owned DAE variables.

## Initialization

From the initialized inputs,

```math
\mathbf{m} \leftarrow \dfrac{2}{v_{\mathrm{dc}}}\mathbf{u}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`m` | [-] | Modulation command | $\mathbf{m} \in \mathbb{R}^3$
