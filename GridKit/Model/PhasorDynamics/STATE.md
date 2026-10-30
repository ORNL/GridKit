# State Input Format

A state file separates the model definition from its initial state.

This also separates the case from the dispatch. One canonical `texas.case.json` can have multiple operating points in `examples/`:

- `texas-summer2025.state.json`
- `texas-winter2025.state.json`
- `texas-summer2026.state.json`

For example `texas-summer2025.state.json` might be structured as:

```json
"header":{
    "format_version": 0.2,
    "format_revision": 0,
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

Buses own their voltage and shunt outputs. Devices own their terminal current outputs. These values support initialization across PowerFlow, OPF, PhasorDynamics, and EMT.

## Bus

### RMS State (PowerFlow & PhasorDynamics)

RMS voltages use real and imaginary parts:

```json
"bus_id_2533":{
    "Vr": 0.9289638822595822,
    "Vi": -0.39534548980249884
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

Device entries give initial outputs and settings such as `in_service`, `tap`, and `phase`.

`in_service` is an optional Boolean operating setting for branches and
injection devices. It defaults to `true`. An out-of-service device
contributes no terminal current to its connected buses.

Use `ir`/`ii` for one terminal and `ir1`/`ii1`, `ir2`/`ii2` for multiple terminals. Numbers follow the model's terminal numbering. Currents are positive into each connected bus. Phasor currents use the system base. EMT uses the model's phase output names and units.

> At nonzero voltage, convert terminal power to initial current with `I = conj((P + jQ) / V)`.

```json
"gen_id_2":{
    "in_service": false,
    "ir": 0.0,
    "ii": 0.0
},
"br_id_2":{
    "in_service": true,
    "tap": 1,
    "phase": 0
}
```

# Migration

- Initialize components from output values and connected inputs.
- Derive missing values and check the initialization equations.
- Each model's `<Model>InitialValues` enum lists the initial outputs and operating settings it accepts.

`Bus`:

- Treat buses as components with `vr` and `vi` output ports.
- Move voltage initialization from the case to the state file.

`Branch`:

- Expose `ir1`, `ii1`, `ir2`, and `ii2` as output ports.
- Make `tap`, `phase`, and `in_service` input ports.

`REGCA`, `LoadZ`, `LoadZIP`, `GENROU`, `GENSAL`, `GenClassical`:

- Expose `ir` and `ii` as output ports.
- Initialize from terminal voltage and current.
- Remove generator and converter `p0`/`q0` parameters.
- Make `in_service` an input port.
