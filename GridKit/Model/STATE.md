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

Buses own their voltage and shunt outputs. Devices own their terminal current outputs. These values support initialization across PowerFlow, OPF, PhasorDynamics, and EMT.

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

Use `ir`/`ii` for one terminal and `ir1`/`ii1`, `ir2`/`ii2` for multiple terminals. Numbers follow the model's terminal numbering. Currents are positive into each connected bus. Phasor currents use the system base. EMT uses the model's phase output names and units.

> At nonzero voltage, convert terminal power to initial current with `I = conj((P + jQ) / V)`.

```json
"gen_id_2":{
    "online": true,
    "ir": 0.8,
    "ii": -0.2
},
"br_id_2":{
    "open": false,
    "tap": 1,
    "phase": 0
}
```

# Migration

- Initialize components from output values and connected inputs.
- Derive missing values and check the initialization equations.

`Bus`:

- Treat buses as components with `vr` and `vi` output ports.
- Move voltage initialization from the case to the state file.

`Branch`:

- Expose `ir1`, `ii1`, `ir2`, and `ii2` as output ports.
- Make `tap`, `phase`, and `open` input ports.

`REGCA`, `LoadZ`, `LoadZIP`, `GENROU`, `GENSAL`, `GenClassical`:

- Expose `ir` and `ii` as output ports.
- Initialize from terminal voltage and current.
- Remove generator and converter `p0`/`q0` parameters.
- Make `online` an input port.

# Readers and Writers

`GridKit/Model/StateData.hpp` reads and writes this format. Device keys are the
case `id`. Bus keys are `bus_id_` followed by the bus number.

- `PhasorDynamics::applyState` writes a state into parsed case data: bus
  `Vr0`/`Vi0`, machine `p0`/`q0`, `LoadZIP` `Pnom`/`Qnom`, and branch `tap`
  and `phase`. Open branches and offline `LoadZ` devices are removed, and an
  offline `LoadZIP` draws no power. Offline machines are rejected.
- `PhasorDynamics::extractState` returns the operating point stored in a case:
  bus voltages, machine and load currents, and branch `tap` and `phase`.
- `OptimalDispatch` reads a state as its starting point and writes the solved
  bus voltages and device currents.
