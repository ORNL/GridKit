# Signal Model

`Signal` holds $N$ signals $\mathbf{s}$ that models read at their inputs. A
held signal keeps the value that the state gives, until an event changes it.
A signal that is not held follows an output: the model that computes the
output adds it to the residual of $\mathbf{s}$. The system adds one `Signal`
for the inputs that a case leaves out, and one for the outputs that inputs
read.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of signals | Required, positive integer
$\mathrm{held}$ | [-] | `held` | Whether the signals are held | Default `false`

### Parameter Validation

```math
N \in \mathbb{Z}_{>0}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{s}$ | `s` | Output | $[s]$ | Signal | Internal variable, $\mathbf{s} \in \mathbb{R}^N$

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{s}$ | $[s]$ | Signal | If held, $\mathbf{s} \in \mathbb{R}^N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{s}$ | $[s]$ | Signal | If not held, $\mathbf{s} \in \mathbb{R}^N$

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

### Internal Equations

#### Differential

If held,

```math
0 = -\dfrac{\mathrm{d}\mathbf{s}}{\mathrm{d}t}
```

#### Algebraic

If not held, with $\mathbf{y}$ the output that the signal follows,

```math
0 = -\mathbf{s} + \mathbf{y}
```

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`s` | $[s]$ | Signal | $\mathbf{s} \in \mathbb{R}^N$
