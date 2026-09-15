# State Input Format

A state file separates the model definition from its initial state.

This also separates the case from the dispatch. One canonical `texas.case.json` can have multiple operating points in `examples/`:

- `texas-summer2025.state.json`
- `texas-winter2025.state.json`
- `texas-summer2026.state.json`

For example `texas-summer2025.state.json` might be structured as:

```json
"header":{
    ...
},
"buses":{
    ...
},
"devices":{
    ...
}
```

Each entry gives initial outputs and operating settings. Models use these values and connected inputs to initialize their own state.

Buses own their voltage and shunt outputs. Devices own their current and power outputs. These values support initialization across PowerFlow, OPF, PhasorDynamics, and EMT.

## Bus

### RMS State (PowerFlow & PhasorDynamics)

RMS voltages use real and imaginary parts:

```json
"bus_id_2533":{
    "vr": 0.9289638822595822,
    "vi": -0.39534548980249884
}
```

### ABC State (EMT)

EMT bus entries give instantaneous phase voltages. Device entries give currents. Each model defines any other initial states or history it needs.

```json
"bus_id_2533":{
    "va": 0.9289638822595822,
    "vb": -0.39534548980249884,
    "vc": -0.39534548980249884
}
```

## Devices

Device entries give initial outputs and settings such as `online`, `open`, `tap`, and `phase`.

Use `ir`/`ii` for one terminal and `ir1`/`ii1`, `ir2`/`ii2` for multiple terminals. Numbers follow the model's terminal numbering. Currents are positive into each connected bus. EMT uses the model's phase output names.

> A model may accept terminal `p`/`q` or `ir`/`ii`. Bus voltage relates the two; if both are supplied, they must agree. Here `p` and `q` mean terminal power. Nominal load values and control references are separate quantities.

```json
"gen_id_2":{
    "online": true,
    "p": 55,
    "q": -12
},
"br_id_2":{
    "open": false,
    "tap": 1,
    "phase": 0
}
```

# Migration

- Initialize each model from its supplied outputs and connected inputs.
- Derive missing values and check supplied values against the model's initialization equations.

`Bus`:

- Remove bus `init` fields from `INPUT_FORMAT.md`.

`Branch`:

- Remove `tap` and `phase` as parameters.
- Add `setTap(..)` and `setPhase(..)`.
- Add `setOpen(..)`.

`LoadZIP`, `GENROU`, `GENSAL`, `GenClassical`:

- Move dispatch values from case parameters to state entries.
- Add `setOnline(..)`.
