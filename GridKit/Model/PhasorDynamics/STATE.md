# State Input Format

A representation of the models state from which it can be initialized. This decouples the model definition from the model state.

In addition to making model initialization easier, it also lets us seperate the case from the dispatch. With one canonical `texas.case.json` we might have multiple operating points provided in `examples/`:
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

Inside the `buses` and `devices` the relevant boundary states needed for intiailziation are given for each component. Then Model initialization at a higher or lower fidelity `PowerFlow`$\leftrightarrow$ `PhasorDynamics`$\leftrightarrow$ `EMT` requires no direct API.

## Bus

### RMS State (PowerFlow & PhasorDynamics)

RMS state initialization/storage is written with real and imaginary parts:

```json
"bus_id_2533":{
    "vr": 0.9289638822595822,
    "vi": -0.39534548980249884,
    "injections": {
        "6533_C":{
            "ir": 13.185009956359863,
            "ii": 0.785739004611969
        },
        "6533_G":{
            "ir": 12.132010459899902,
            "ii": 0.785739004611969
        }
    }
},
```

### ABC State (EMT)

EMT should support balanced initialization (see previous section) but also unbalanced initialization, which would require higher fidelity. This would also specifiy branch injections for distirbuted network models.

```json
"bus_id_2533":{
    "va": 0.9289638822595822,
    "vb": -0.39534548980249884,
    "vc": -0.39534548980249884,
    "injections": {
        "gen_325":{
            "ia": 13.185009956359863,
            "ib": 0.785739004611969,
            "ic": 0.785739004611969
        },
        "branch_6533":{
            "ia": 13.185009956359863,
            "ib": 0.785739004611969,
            "ic": 0.785739004611969
        }
    }
},
```

## Devices

This is where we store stateful information like connectivity and dispatch, which can be shared across all domains (PowerFlow & PhasorDynamics & EMT). This would also let us specify non steady state initialziations with `"omega": 0.0`


```json
"gen_id_2":{
    "online": true,
    "p": 55,
    "q": -12
}
"br_id_2":{
    "open": false,
    "tap": 1,
    "phase": 0
}
```


# Migration

`Bus`:
- Remove the `init` fields from `Bus` in the `INPUT_FORMAT.md` 
- make an interface `setVr(..)` and `setVi(..)`.

`Branch`:
- Remove `tap` and `phase` as parameters
- make an interface `setTap(..)` and `setPhase(..)`.
- add `setOpen(..)` interface 

`LoadZIP`, `GENROU`, `GENSAL`, `GenClassical`:
- Remove `p0` and `q0` as parameters
- make an interface `setP(..)` and `setQ(..)`.
- add `setOnline(..)` interface 
