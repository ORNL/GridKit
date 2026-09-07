# Signal

Signals provide scalar connection points between EMT components.
Components attach external inputs to signals and assign internal outputs
to signals through `ComponentSignals`.

A linked signal stores a pointer to the component variable that owns the
signal value, so other connected components can read or initialize that value
without owning the producing model. A fully bound signal additionally
stores pointers to the owning variable's derivative and residual row, so
connected components can read the derivative and accumulate external residual
contributions into the owner's residual row. Three-phase connections use
three scalar signals, one for each phase.

Computed algebraic signals bind a value getter and its gradient with respect
to global DAE variables. They own no variable index, derivative, or residual
row, and cannot be initialized by a consumer. Reads evaluate the current
inputs recursively. Residual Jacobians compose the local input derivatives
with the signal gradients. The input graph must be acyclic, and its gradient
structure must remain fixed after allocation, including entries whose current
coefficient is zero.

A declared constant (`{"id": "pref", "value": 0.2}`) owns its finite value
without a DAE variable, derivative, residual row, or gradient entries.
`bindConstant(value)` reserves its producer; `setConstantValue(value)` updates
only a declared constant. Consumer initialization cannot overwrite it.
Rebinding a signal to a component variable or expression clears its constant
status. The EMT application's `signal_step` event uses this explicit update
operation to change references at a scheduled time.

## Model Parameters

Symbol | Description
-------|------------
`id` | Unique string identifier for the signal
`value` | Optional initial value of a declared, externally writable constant
