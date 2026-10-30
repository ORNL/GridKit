# GasBus

A gas bus owns the pressure $p$ and the mass-flow balance residual. Each
connected device adds its mass flow after the bus resets the residual to the
injection.

## Notes

- Mass flow entering the bus has positive sign, so the injection $q$ of a gas
  load is negative.
- Connected GasBranch models add their terminal storage, so $p$ is a
  differential variable.
- A `withdrawal` event replaces $q$ with the negative of the withdrawal, for
  example for a leak.

## Model Parameters

Symbol | Units  | JSON | Description | Typical Value | Note
-------|--------|------|-------------|---------------|-----
$q$    | [kg/s] | `q`  | Injection   | 0.0           |

### Parameter Validation

```math
|q| < \infty
```

### Model Derived Parameters

None.

## Model Ports

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description  | Note
-------|-------|--------------|-----
$p$    | [Pa]  | Bus pressure |

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

Let $\mathcal{D}$ denote the set of components connected to the bus, and let
$q_d$ be the mass flow component $d$ adds.

### Internal Equations

#### Differential

```math
0 = q + \sum_{d \in \mathcal{D}} q_d
```

#### Algebraic

None.

### External Equations

None.

## Initialization

### Internal Initialization

The initial pressure must be positive.

```math
p \leftarrow \text{initial bus pressure}
```

The derivative initializes to zero.

## Monitors

Monitor | Units | Description  | Note
--------|-------|--------------|-----
`p`     | [Pa]  | Bus pressure |
