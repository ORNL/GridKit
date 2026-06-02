# Input file for GridKit phasor dynamics application

`DynamicSimulation` accepts event studies and forced-oscillation studies. Event
studies use `events`. Forced-oscillation studies use `forced_oscillations` and
must not include `events`.

## Root elements

   Name               | Value
 ---------------------|-------------------------------------------------------
  `system_model_file` | Path to the system model file[^1]
  `dt`                | A floating point value for time step size
  `tmax`              | A floating point value for max time
  `events`            | Event-study array of event groups (see [Events](#events) below)
  `forced_oscillations` | Forced-oscillation-study injection array (see [Forced Oscillations](#forced-oscillations) below)
  `max_step`          | Optional IDA internal max-step cap
  `output_file`       | Path to output (CSV) file
  `reference_file`    | A string containing the name of the case
  `error_tolerance`   | A string containing the name of the case

[^1]: See system model [case format](../../Model/PhasorDynamics/INPUT_FORMAT.md)

## Events

Each event group describes a system event that occurs at a given time point

   Name              | Value
 --------------------|-------------------------------------------------------
  `time`             | A floating point value for time event occurs
  `type`             | Event type (one of { "fault_on", "fault_off" })
  `element_id`       | An integer value referencing the element associated with the event (e.g., bus fault id)

## Forced Oscillations

Each forced-oscillation entry creates a `ForcedOscillation` signal operator in
the parsed model before `SystemModel` is constructed. Injections perturb the DAE.

   Name     | Value
 -----------|-------------------------------------------------------
  `id`      | Unique forced-oscillation device ID
  `target`  | `Class:id.port` or `signal:ID`
  `mode`    | One of { `add`, `drive`, `signal` }
  `params`  | `ForcedOscillation` parameters such as `A`, `f`, `Kf`, `Phi`, `Ton`
  `mon`     | Optional monitor variables; defaults to `in`, `force`, `out`, `active`

Modes:

- `add`: the target must be an existing consumer signal port. The old signal is
  the forced-oscillation input, and the target port is retargeted to the forced
  output.
- `drive`: the target must be a consumer signal port. The forced oscillation
  drives that port directly, including optional ports absent from the base case.
- `signal`: the target is an explicit signal ID. The source may be floating if
  it is only being monitored.

`max_step` should be small enough to resolve the injected waveform. A practical
starting point is about 20 samples per forcing cycle:

```text
dt <= 1 / (20 * f_max)
```

When `max_step` is omitted in a forced-oscillation study, the application uses
`min(dt, 1 / (20 * f_max))`; if `f_max` is zero, it uses `dt`.

Example:

```json
{
  "system_model_file": "TwoBusTgov1.json",
  "dt": 0.00416666666666,
  "tmax": 20.0,
  "output_file": "TwoBusTgov1_fo.csv",
  "forced_oscillations": [
    {
      "id": "fo_pmech",
      "target": "Genrou:DV1.pmech",
      "mode": "add",
      "params": { "A": 0.02, "f": 0.7, "Ton": 2.0, "Toff": 15.0, "Tr": 0.5, "Tf": 0.5 },
      "mon": ["in", "force", "out", "active"]
    }
  ]
}
```
