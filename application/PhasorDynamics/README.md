# Input file for GridKit phasor dynamics application

## Root elements

   Name               | Value
 ---------------------|-------------------------------------------------------
  `system_model_file` | Path to the system model file[^1]
  `dt_monitor`        | Monitor output time interval for recorded simulation results (default: 0, no intermediate monitoring)
  `tmax`              | A floating-point value for max time
  `rel_tol`           | Relative solver tolerance (default: 1.0e-7)
  `abs_tol`           | Absolute solver tolerance override (default: 1.0e-9)
  `dt_fixed`          | Fixed solver time step size, or 0 for adaptive stepping (default: 0)
  `max_steps`         | Maximum number of solver time steps, 0 for the IDA default, or a negative number for unlimited steps (default: 0)
  `events`            | An array of event groups (see [Events](#events) below)
  `forced_oscillations` | Optional array of additive forced-oscillation injections (see [Forced Oscillations](#forced-oscillations) below)
  `output_file`       | Path to output (CSV) file (optional)
  `reference_file`    | A string containing the name of the case (optional)
  `error_type`        | One of { "relative" (default), "absolute" }
  `error_tolerance`   | A floating-point value for highest allowable total error (default: 1.0e-4)
  `abs_err_threshold` | A floating-point value for the smallest value at which to scale relative error (default: machine epsilon for double-precision)

[^1]: See system model [case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)

## Forced Oscillations

For copy-ready study files and common waveform recipes, see the
[forced-oscillation study guide](FO.md).

Each `forced_oscillations` entry adds a `ForcedOscillation` signal source and
an algebraic signal-node junction before `SystemModel` is constructed. The
source waveform $s_{\mathrm{FO}}$ is added to the target's existing input $u$,

```math
y = u+s_{\mathrm{FO}},
```

and the target input is redirected to junction output $y$. The original input
is the junction's `initialization_input`, so initialization writes from the
target are propagated back through the pre-existing signal path. Study
injections and `events` may be used together; an empty `events` array isolates
the forced response. See the
[signal-source documentation](../../GridKit/Model/PhasorDynamics/SignalSource/README.md)
for the complete waveform equation.

### Injection Fields

| Name | Description | Default |
| --- | --- | --- |
| `id` | Unique `ForcedOscillation` component ID | Required |
| `target` | Existing consumer input in `Class:id.port` form | Required |
| `mode` | Injection mode; only `add` is supported | `add` |
| `params` | Object containing source parameters | Empty object |
| `mon` | Source variables to record | `output`, `envelope`, `active` |

`Class`, `id`, and `port` are case-sensitive. `port` must name an entry in the
target model's signal-input set, and that input must already reference a signal.
Output or electrical-bus ports cannot be targeted. Monitor names must be one of
`output`, `envelope`, or `active`. Injection IDs must be nonempty and unique;
each target may occur at most once in a study.

### Source Parameters

The target signal defines the units and per-unit base, if any, of the source
output and amplitude $A$.

| Symbol | Units | JSON | Description | Default |
| --- | --- | --- | --- | --- |
| $A$ | Target-signal units | `A` | Oscillation amplitude | 0.0 |
| $f$ | Hz | `f` | Initial oscillation frequency | 0.0 |
| $K_f$ | Hz/s | `Kf` | Linear frequency ramp | 0.0 |
| $\Phi$ | rad | `Phi` | Phase offset | 0.0 |
| $T_{\mathrm{on}}$ | s | `Ton` | Activation time | 0.0 |
| $T_{\mathrm{off}}$ | s | `Toff` | Deactivation time; a negative value disables deactivation | -1.0 |
| $T_r$ | s | `Tr` | Raised-cosine rise time | 0.0 |
| $T_f$ | s | `Tf` | Raised-cosine fall time | 0.0 |
| $K_d$ | 1/s | `Kd` | Exponential decay rate | 0.0 |
| $W$ | integer | `waveform` | Carrier selector: 0 sine, 1 square, 2 triangle, or 3 sawtooth | 0 |

All supplied values must be finite. $A$, $f$, $K_f$, $T_r$, $T_f$, and $K_d$
must be non-negative. A non-negative $T_{\mathrm{off}}$ must satisfy
$T_{\mathrm{off}} \ge T_{\mathrm{on}}$. The waveform selector must be an
exact integer from 0 through 3.

### Monitorable Outputs

| Name | Description |
| --- | --- |
| `output` | Published waveform $s_{\mathrm{FO}}$ |
| `envelope` | Raised-cosine activation envelope |
| `active` | Active-window indicator, 0.0 or 1.0 |

### Example

```json
{
  "system_model_file": "wecc.json",
  "dt_monitor": 0.00416666666666,
  "tmax": 30.0,
  "events": [],
  "forced_oscillations": [
    {
      "id": "fo_1032_C_pmech",
      "target": "Genrou:1032_C_genrou.pmech",
      "mode": "add",
      "params": {
        "A": 0.1,
        "f": 0.7,
        "Ton": 10.0,
        "Toff": -1.0,
        "Tr": 1.0
      },
      "mon": ["output", "envelope", "active"]
    }
  ]
}
```
