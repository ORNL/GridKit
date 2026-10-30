# GasBus

A gas bus owns the pressure $p$ and the mass-flow balance residual. Each
connected device adds its mass flow after the bus resets the residual to zero.

## Notes

- Mass flow entering the bus has positive sign.
- Connected GasBranch models add their terminal storage, so $p$ is a
  differential variable.

## Model Parameters

None.

### Parameter Validation

None.

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
0 = \sum_{d \in \mathcal{D}} q_d
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
