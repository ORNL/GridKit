# Forced-oscillation studies

`ForcedOscillationStudy.hpp` adds forced-oscillation injections to the study
file consumed by `DynamicSimulation`; there is no separate study executable.
Each target receives `original_input + forcing` without changing the case file
on disk.

## Minimal study

For example, save this as
`examples/PhasorDynamics/validation/WECC240/my-fo.solver.json`:

```json
{
  "system_model_file": "../../../../cases/PhasorDynamics/WECC240/WECC240.case.json",
  "dt_monitor": 0.05,
  "tmax": 10.0,
  "events": [],
  "forced_oscillations": [
    {
      "id": "fo_1032_C_pmech",
      "target": "Genrou:1032_C_genrou.pmech",
      "mode": "add",
      "params": {
        "A": 0.02,
        "f": 0.7,
        "Ton": 1.0,
        "Toff": -1.0,
        "Tr": 0.15,
        "waveform": 0
      },
      "mon": ["output", "envelope", "active"]
    }
  ]
}
```

Build and run from the repository root:

```sh
cmake --build build --target DynamicSimulation -j 10
build/application/PhasorDynamics/DynamicSimulation \
  examples/PhasorDynamics/validation/WECC240/my-fo.solver.json
```

`system_model_file` is resolved relative to the `.solver.json` file.
`dt_monitor` controls output sampling, not the solver time step. Keep
`events` as an empty array to isolate the forced response; ordinary events and
forced oscillations may also be used in the same study.

## Injection settings

This object shows every injection field and every source parameter:

```json
{
  "id": "unique_source_id",
  "target": "Class:component_id.signal_input",
  "mode": "add",
  "params": {
    "A": 0.02,
    "f": 0.7,
    "Kf": 0.0,
    "Phi": 0.0,
    "Ton": 1.0,
    "Toff": -1.0,
    "Tr": 0.15,
    "Tf": 0.0,
    "Kd": 0.0,
    "waveform": 0
  },
  "mon": ["output", "envelope", "active"]
}
```

| Key | Meaning | Default |
| --- | --- | --- |
| `A` | Peak amplitude in the target signal's units and per-unit base | `0.0` |
| `f` | Initial frequency in Hz | `0.0` |
| `Kf` | Linear frequency ramp in Hz/s | `0.0` |
| `Phi` | Phase offset in radians | `0.0` |
| `Ton` | Activation time in seconds | `0.0` |
| `Toff` | Deactivation time; a negative value keeps the source on | `-1.0` |
| `Tr` | Raised-cosine rise time in seconds | `0.0` |
| `Tf` | Raised-cosine fall time in seconds | `0.0` |
| `Kd` | Exponential decay rate in 1/s | `0.0` |
| `waveform` | `0` sine, `1` smooth square, `2` smooth triangle, `3` smooth sawtooth | `0` |

All supplied values must be finite JSON numbers. `A`, `f`, `Kf`, `Tr`, `Tf`,
and `Kd` must be nonnegative. `waveform` must be an exact integer from 0
through 3. A nonnegative `Toff` must be greater than or equal to `Ton`; `Tf`
is used only with a nonnegative `Toff`.

### Common waveforms

Replace an injection's `params` object with any of these examples.

Stationary sine that remains active through the simulation:

```json
{"A": 0.02, "f": 0.7, "Ton": 1.0, "Toff": -1.0, "Tr": 0.15, "waveform": 0}
```

Nonstationary chirp; its instantaneous frequency after activation is
`f + Kf * (t - Ton)`:

```json
{"A": 0.02, "f": 0.2, "Kf": 0.15, "Ton": 1.0, "Toff": -1.0, "Tr": 0.15, "waveform": 0}
```

Exponentially decaying oscillation:

```json
{"A": 0.04, "f": 0.7, "Kd": 0.25, "Ton": 1.0, "Toff": -1.0, "Tr": 0.15, "waveform": 0}
```

Finite smooth square-wave injection:

```json
{"A": 0.02, "f": 0.7, "Phi": 0.0, "Ton": 1.0, "Toff": 8.0, "Tr": 0.15, "Tf": 0.25, "waveform": 1}
```

Use the same settings with `"waveform": 2` for a smooth triangle or
`"waveform": 3` for a smooth sawtooth.

## Targets and multiple injections

`target` is case-sensitive and has the form `Class:id.port`. The port must be
an existing signal input in the target device's `ports` object. Electrical-bus
ports and signal outputs cannot be targeted. For example, this array injects
at one synchronous-machine control and one renewable-converter control:

```json
{
  "forced_oscillations": [
    {
      "id": "fo_sync_pmech",
      "target": "Genrou:1032_C_genrou.pmech",
      "params": {"A": 0.02, "f": 0.7, "Ton": 1.0, "Toff": -1.0, "waveform": 0},
      "mon": ["output"]
    },
    {
      "id": "fo_renewable_ipcmd",
      "target": "Regca:1032_S_regca.ipcmd",
      "params": {"A": 0.01, "f": 0.4, "Kf": 0.1, "Ton": 1.0, "Toff": -1.0, "waveform": 2},
      "mon": ["output", "envelope", "active"]
    }
  ]
}
```

Injection IDs and targets must each be unique within a study. `mode` may be
omitted because it defaults to `add`; no other mode is supported. `mon` may be
omitted to record all three source diagnostics, or set to an empty array to
record none.

## Recording the system response

The injection's `mon` field records only the source. Add response monitors to
the existing objects in the case file. For example, merge the shown `mon`
fields into the complete bus and device objects:

```json
{
  "monitors": [{"file_name": "response.csv", "format": "CSV"}],
  "buses": [
    {"number": 1032, "mon": ["Vm", "Va"]}
  ],
  "devices": [
    {"class": "Genrou", "id": "1032_C_genrou", "mon": ["p", "q", "omega"]},
    {"class": "Regca", "id": "1032_S_regca", "mon": ["p", "q"]}
  ]
}
```

The case file's `monitors` array normally owns the output sink. The optional
solver-file `output_file` is most useful when the case has no CSV sink; an
existing case CSV sink takes precedence.

## Optional solver controls

These keys belong at the root of the `.solver.json` file:

```json
{
  "rel_tol": 1e-7,
  "abs_tol": 1e-9,
  "dt_fixed": 0.0,
  "max_steps": 0,
  "output_file": "response.csv",
  "reference_file": "reference.csv",
  "error_type": "relative",
  "error_tolerance": 1e-4,
  "abs_err_threshold": 1e-12
}
```

`dt_fixed: 0.0` selects adaptive stepping; a positive value selects fixed
stepping. `max_steps: 0` uses the IDA default. The reference and error fields
are needed only when comparing the run against reference data;
`error_tolerance` may be a number or an array. Relative `reference_file` paths
are resolved from the solver-file directory, while relative `output_file`
paths are resolved from the process working directory.

See the complete [WECC240 example][wecc-example],
the [waveform equations](../../GridKit/Model/PhasorDynamics/SignalSource/README.md),
and the [100-study gallery generator](../../examples/FO/generate_gallery.py).

[wecc-example]: ../../examples/PhasorDynamics/validation/WECC240/WECC240.forced-oscillation.solver.json
