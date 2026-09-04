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

## Model Parameters

Symbol | Description
-------|------------
`signal_id` | Unique identifier for the signal
