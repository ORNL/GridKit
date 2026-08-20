# Signal Node

Signal nodes provide scalar connection points between phasor-dynamics
components. Components attach external inputs to signal nodes and assign
internal outputs to signal nodes through `SignalPorts`.

A linked signal node stores a pointer to the component variable that owns the
signal value, so other connected components can read or initialize that value
without owning the producing model.

## Model Parameters

Symbol | Description
-------|------------
`signal_id` | Unique identifier for the signal node
