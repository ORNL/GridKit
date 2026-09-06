# SignalNode

Signal nodes provide scalar connection points between phasor-dynamics
components. Components attach external inputs to signal nodes and assign
internal outputs to signal nodes through `ComponentSignals`.

A linked signal node stores a pointer to the component variable that owns the
signal value, so other connected components can read or initialize that value
without owning the producing model.

## Model Parameters

Symbol      | Units | JSON        | Description                           | Typical Value | Note
------------|-------|-------------|---------------------------------------|---------------|-----
`signal_id` | [-]   | `signal_id` | Unique identifier for the signal node |               |

### Parameter Validation

None.

### Model Derived Parameters

None.

## Model Ports

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

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

None.

## Initialization

None.

## Monitors

None.
