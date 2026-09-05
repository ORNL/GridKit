# Signal

Signals provide scalar connection points between EMT components.
Components attach external inputs to signals and assign internal outputs
to signals through `ComponentSignals`.

A linked signal stores a pointer to the component variable that owns the
signal value, so other connected components can read or initialize that value
without owning the producing model. A fully bound signal additionally
stores pointers to the owning variable's derivative and residual row, so
connected components can read the derivative and accumulate external residual
contributions into the owner's residual row. Three fully bound signals
form a `Port3`, the three-phase electrical connection point owned by the
component that owns the phase variables.

Computed algebraic signals bind a value getter and its gradient with respect
to global DAE variables. They own no variable index, derivative, or residual
row, and cannot be initialized by a consumer. Reads evaluate the current
inputs recursively. Residual Jacobians compose the local input derivatives
with the signal gradients. The input graph must be acyclic, and its gradient
structure must remain fixed after allocation, including entries whose current
coefficient is zero.

## Model Parameters

Symbol | Description
-------|------------
`id` | Unique string identifier for the signal
