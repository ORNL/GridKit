# Propagation Model

`Propagation` represents an EMT modal propagation function.

```math
\mathbf{H}(s) =
\sum_{m=1}^{M} e^{-s\tau_m}\mathbf{H}^{\mathrm{mps}}_m(s)
```

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../docs/Figures/EMT/Propagation/diagram.png">

  Figure 1: Propagation model
</div>

## Model Parameters

Symbol                | Units | JSON  | Description              | Typical Value | Note
----------------------|-------|-------|--------------------------|---------------|-----
$\boldsymbol{\tau}$   | s     | `tau` | Modal propagation delays | --            | $\boldsymbol{\tau} \in \mathbb{R}^M$

### Parameter Validation

```math
\tau_m > 0,\qquad m = 1,\dots,M
```

The number of modal `StateSpace` submodels, `Delay` submodels, and entries in
`tau` must match. Each submodel must satisfy its own validation rules. All
`StateSpace` submodels must accept the same input dimension and produce the same
output dimension.

### Model Derived Parameters

None.

## Model Submodels

Submodel     | Qty  | Symbol                              | Description                         | Note
-------------|------|-------------------------------------|-------------------------------------|-----
`StateSpace` | $M$  | $\mathbf{H}^{\mathrm{mps}}_m(s)$    | Modal propagation shaping function  | $\mathbf{H}^{\mathrm{mps}}_m(s) \in \mathbb{R}^{K \times K}$
`Delay`      | $KM$ | $e^{-s\tau_m}$                      | Modal propagation delay             |

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol       | Units | Description              | Note
-------------|-------|--------------------------|-----
$\mathbf{y}$ | [-]   | Propagation output state | $\mathbf{y} \in \mathbb{R}^K$

### External Variables

#### Differential

Symbol         | Units | Description                         | Note
---------------|-------|-------------------------------------|-----
$\mathbf{u}$   | [-]   | Input vector read from `input` port | $\mathbf{u} \in \mathbb{R}^K$

#### Algebraic

None.

## Model Equations

### Differential Equations

None.

### Algebraic Equations

With $\mathbf{w}_m$ the output of the $m$th modal `Delay` submodel:

```math
0 = -\mathbf{y} + \sum_{m=1}^{M}\mathbf{w}_m
```

## Initialization

Submodels initialize according to their own model contracts. The propagation
output initializes from the modal `Delay` outputs:

```math
\mathbf{y}_0 = \sum_{m=1}^{M}\mathbf{w}_{m,0}
```

## Model Outputs

Output | Units | Description              | Note
-------|-------|--------------------------|-----
`out`  | [-]   | Propagated output vector | $\mathbf{y} \in \mathbb{R}^K$
