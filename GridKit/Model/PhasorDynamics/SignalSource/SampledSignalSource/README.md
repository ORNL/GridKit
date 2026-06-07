# **Sampled Signal Source Model**

The Sampled Signal Source model writes an exogenous sampled waveform to a
signal node. The model has no solver-owned state or residual equations.

Notes:
- Samples are linearly interpolated.
- Values outside the sampled time range use the nearest endpoint value.
- CSV file paths are resolved relative to the system model JSON file.

## Model Parameters

Symbol   | Units | JSON     | Description                             | Typical Value
---------|-------|----------|-----------------------------------------|--------------
$k$      | [-]   | `scale`  | Multiplicative factor on sampled values | 1
$b$      | [-]   | `offset` | Additive offset after scaling           | 0

### Source Data

The source data is provided through the device `source` object.

Field          | Description
---------------|------------
`type`         | `samples` for inline data or `csv` for file input
`samples`      | Inline array of `[time, value]` pairs
`file`         | CSV file path
`time_column`  | CSV column containing sample times
`value_column` | CSV column containing sample values

### Parameter Validation

Invalid source data is rejected by the following checks:

```math
\begin{aligned}
  &\text{type} \in \{\text{samples}, \text{csv}\} \\
  &t_i < t_{i+1}
\end{aligned}
```

At least one sample is required.

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

Symbol | Units | Description
-------|-------|------------
$u_s$  | [-]   | Signal value written to the `output` port

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

## Initialization

The output signal initializes from the sampled source value at the initial
simulation time.

## Model Outputs

Output | Units | Description
-------|-------|------------
`out`  | [-]   | Sampled source signal

For sample interpolation $s(t)$, the output signal is:

```math
u_s(t) = k s(t) + b
```
