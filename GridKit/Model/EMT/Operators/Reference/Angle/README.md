# Angle Model

`Angle` integrates electrical angular frequency to obtain the reference angle
used by the [Park transformation](../Park/README.md). The angle is continuous
and is not wrapped at multiples of $2\pi$.

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
$\omega$ | `omega` | Input | [rad/s] | Electrical angular frequency | Required
$\theta$ | `theta` | Output | [rad] | Electrical reference angle | d-axis relative to the phase-a axis

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\theta$ | [rad] | Electrical reference angle |

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\omega$ | [rad/s] | Electrical angular frequency | Differential-input configuration

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\omega$ | [rad/s] | Electrical angular frequency | Algebraic-input configuration

## Model Equations

### Internal Equations

#### Differential

```math
0 = -\dfrac{\mathrm{d}\theta}{\mathrm{d}t} + \omega
```

#### Algebraic

None.

### External Equations

None.

## Initialization

The state-file key `theta` supplies a finite angle, defaulting to zero. The
consistent-initial-condition solve preserves this angle and obtains its
derivative from the connected frequency:

```math
\dfrac{\mathrm{d}\theta}{\mathrm{d}t} \leftarrow \omega
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`theta` | [rad] | Electrical reference angle |
