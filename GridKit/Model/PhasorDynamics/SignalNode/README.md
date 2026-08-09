# Signal Node

Signal nodes provide scalar connection points between phasor-dynamics
components. Components attach external inputs to signal nodes and assign
internal outputs to signal nodes through `ComponentSignals`.

## Ordinary Signal Nodes

A linked signal node stores a pointer to the component variable that owns the
signal value, so other connected components can read or initialize that value
without owning the producing model.

## Algebraic Junctions

A signal node may instead define a weighted algebraic junction. For output
signal $y$, bias $b$, input signals $u_i$, and gains $g_i$, the junction adds
one algebraic equation

```math
f = y - b - \sum_i g_i u_i = 0.
```

The junction owns $y$ as an algebraic DAE variable. Its inputs may be ordinary
signals, exogenous source signals, or other junction outputs. Junctions are
ordered by their dependencies; direct or indirect junction cycles are invalid.

### Initialization

`initialization_input` is the signal ID of the input that receives backward
initialization writes. If a consumer initializes the junction output to
$y_0$, the selected input $u_k$ is initialized as

```math
u_k = \frac{y_0 - b - \sum_{i\ne k} g_i u_i}{g_k}.
```

The selected signal must occur exactly once in `inputs`, and its gain $g_k$
must be nonzero. This designation is a signal ID, not an array position.
The junction output cannot also be assigned by a device signal-output port.

## Model Data

| Name | Description | Default |
| --- | --- | --- |
| `name` | Signal name or description | Required |
| `signal_id` | Unique non-negative signal identifier | Required |
| `junction` | Optional algebraic-junction object | Ordinary signal node |

A `junction` object has the following fields:

| Name | Description | Default |
| --- | --- | --- |
| `bias` | Constant $b$ in output-signal units | 0.0 |
| `initialization_input` | Signal ID selected for backward initialization | Required |
| `inputs` | Nonempty array of weighted input objects | Required |

Each weighted input object has the following fields:

| Name | Description | Default |
| --- | --- | --- |
| `signal_id` | Input signal identifier | Required |
| `gain` | Multiplier $g_i$ applied to the input | 1.0 |

Input signal IDs must be unique within a junction and cannot equal the
junction output's signal ID. The bias and every input gain must be finite.
