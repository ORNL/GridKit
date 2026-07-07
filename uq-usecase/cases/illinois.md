# Illinois case — MATPOWER .m file

## Geo plot

Interactive Plotly map — branch loading (viridis), buses sized by PD, generators colored by fuel type.  
Generated from `m_viz.ipynb` with `show_loading=True`, `fault_bus=None` (parametric — all 200 buses are wired with `BusFault` devices).

[Open map in new tab (larger view) ↗](ACTIVSg200_geo_loading.html){target="_blank"}

<iframe src="ACTIVSg200_geo_loading.html" width="100%" height="650px" style="border:none;"></iframe>

## Gen table (49 generators)

`gen_row` is the 1-based row index in the MATPOWER `mpc.gen` matrix. `PG`/`PMAX` are in MW. Rows marked `*` have `GEN_STATUS = 0` (offline). Fuel types: coal (25), ng (17), wind (6), nuclear (1).

| gen_row | bus | bus_name          | PG (MW) | PMAX (MW) | status | fuel    |
|--------:|----:|:------------------|--------:|----------:|:------:|:--------|
|   1 |  49 | RANTOUL 2 1       |    1.36 |      4.53 | ✓      | coal    |
|   2 |  50 | RANTOUL 2 2       |    1.36 |      4.53 | ✓      | coal    |
|   3 |  51 | RANTOUL 2 3       |    1.36 |      4.53 | ✓      | coal    |
|   4 |  52 | RANTOUL 2 4       |    1.36 |      4.53 | ✓      | coal    |
|   5 |  53 | RANTOUL 2 5       |    2.72 |      9.07 | ✓      | coal    |
|   6 |  65 | PAXTON 1 1        |   86.50 |     86.50 | ✓      | wind    |
|   7 |  67 | MOUNT ZION 1      |    1.41 |      4.70 | ✓      | coal    |
|   8 |  68 | MOUNT ZION 2      |    8.38 |     27.92 | ✓      | coal    |
|   9 |  69 | MOUNT ZION 3      |    8.38 |     27.92 | ✓      | coal    |
|  10 |  70 | MOUNT ZION 4      |    8.38 |     27.92 | ✓      | coal    |
|  11 |  71 | MOUNT ZION 5      |    8.38 |     27.92 | ✓      | coal    |
|  12 |  72 | MOUNT ZION 6      |    8.38 |     27.92 | ✓      | coal    |
|  13 |  73 | MOUNT ZION 7      |    8.38 |     27.92 | ✓      | coal    |
|  14 |  76 | BRIMFIELD 1       |    1.20 |      4.00 | ✓      | ng      |
|  15 |  77 | BRIMFIELD 2       |    0.72 |      2.40 | ✓      | ng      |
|  16 |  78 | BRIMFIELD 3       |    0.00 |     18.00 | ✗ off  | ng      |
|  17 |  79 | BRIMFIELD 4       |    0.00 |     18.00 | ✗ off  | ng      |
|  18 |  90 | CLINTON 3 1       |    0.96 |      3.20 | ✓      | ng      |
|  19 |  91 | CLINTON 3 2       |    1.50 |      5.00 | ✓      | ng      |
|  20 |  92 | CLINTON 3 3       |    0.00 |      6.30 | ✗ off  | ng      |
|  21 |  94 | TUSCOLA 2 1       |    5.40 |     18.00 | ✓      | coal    |
|  22 | 104 | ELLSWORTH 1 2     |   67.60 |     67.60 | ✓      | wind    |
|  23 | 105 | ELLSWORTH 1 3     |  154.80 |    154.80 | ✓      | wind    |
|  24 | 114 | NORMAL 2 2        |    1.40 |      1.40 | ✓      | wind    |
|  25 | 115 | NORMAL 2 3        |  133.50 |    133.50 | ✓      | wind    |
|  26 | 125 | BARTONVILLE 2     |   39.02 |    130.05 | ✓      | coal    |
|  27 | 126 | BARTONVILLE 3     |   39.02 |    130.05 | ✓      | coal    |
|  28 | 127 | BARTONVILLE 4     |   39.02 |    130.05 | ✓      | coal    |
|  29 | 135 | PEKIN 1 2         |  133.92 |    446.40 | ✓      | coal    |
|  30 | 136 | PEKIN 1 3         |  133.92 |    446.40 | ✓      | coal    |
|  31 | 147 | HOPEDALE 2 1      |   92.40 |     92.40 | ✓      | wind    |
|  32 | 151 | SPRINGFIELD 5 2   |    1.62 |      5.40 | ✓      | coal    |
|  33 | 152 | SPRINGFIELD 5 3   |   23.17 |     77.22 | ✓      | coal    |
|  34 | 153 | SPRINGFIELD 5 4   |   23.17 |     77.22 | ✓      | coal    |
|  35 | 154 | SPRINGFIELD 5 5   |   23.17 |     77.22 | ✓      | coal    |
|  36 | 155 | SPRINGFIELD 5 6   |   23.17 |     77.22 | ✓      | coal    |
|  37 | 161 | SPRINGFIELD 4 1   |    0.00 |    138.60 | ✗ off† | ng      |
|  38 | 164 | CHAMPAIGN 1 1     |    0.00 |     12.00 | ✗ off  | ng      |
|  39 | 165 | CHAMPAIGN 1 2     |    0.00 |     26.00 | ✗ off  | ng      |
|  40 | 166 | CHAMPAIGN 1 3     |    0.00 |      9.40 | ✗ off  | ng      |
|  41 | 167 | CHAMPAIGN 1 4     |    2.82 |      9.40 | ✓      | ng      |
|  42 | 168 | CHAMPAIGN 1 5     |    0.00 |      9.40 | ✗ off  | ng      |
|  43 | 169 | CHAMPAIGN 1 6     |    0.00 |      9.40 | ✗ off  | ng      |
|  44 | 170 | CHAMPAIGN 1 7     |    2.82 |      9.40 | ✓      | ng      |
|  45 | 182 | SPRINGFIELD 2 1   |    5.25 |     17.50 | ✓      | coal    |
|  46 | 183 | SPRINGFIELD 2 2   |    7.98 |     26.60 | ✓      | coal    |
|  47 | 189 | CLINTON 1 2       |  384.37 |    569.15 | ✓      | nuclear |
|  48 | 196 | GIBSON CITY 1 1   |    0.00 |     67.50 | ✗ off  | ng      |
|  49 | 197 | GIBSON CITY 1 2   |    0.00 |     67.50 | ✗ off† | ng      |

† Bus 161 and bus 197 have `GEN_STATUS=0` in the `.m` file but **are included** in `illinois.json` as Genrou devices. See [Discrepancy vs. .m file](#discrepancy-vs-m-file) below.

---

# Illinois case illinois.json file
> examples/PhasorDynamics/Large/Illinois/illinois.json

## Buses (200 buses)

200 buses numbered 1–200, with names and base voltages matching the `.m` file's `mpc.bus_name` and `BASE_KV` columns.

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"number"` | `BUS_I` (col 1) | ✓ 1–200, identical |
| `"name"` | `mpc.bus_name` | ✓ identical |
| `"v_base"` | `BASE_KV` (col 10) | ✓ identical |
| `"init": {Vr, Vi}` | `VM`/`VA` cols (powerflow solution) | ✓ matches powerflow initial conditions |
| `"mon": ["Vr", "Vi"]` | — (GridKit-specific monitoring config) | n/a |
| *(absent)* | `BUS_TYPE` (col 2, PQ/PV/ref) | ✗ not in JSON |
| *(absent)* | `PD`/`QD` (cols 3–4, load) | ✗ not in JSON — loads stored as `"class": "LoadZIP"` devices |
| *(absent)* | `GS`/`BS` (cols 5–6, shunt) | ✗ not in JSON |
| *(absent)* | `AREA`, `ZONE` (cols 7, 11) | ✗ not in JSON |
| *(absent)* | `VMAX`/`VMIN` (cols 12–13) | ✗ not in JSON |

## Branches (246 in JSON vs 245 in .m)

246 `Branch` devices in the `"devices"` array. The `.m` file has 245 branches — the JSON contains one additional branch, bus 82 → bus 64, which has no corresponding row in `mpc.branch`. No parallel circuits: all 246 branch pairs are unique from-to combinations.

Branch IDs follow the pattern `BR_<from>_<to>_<parallel_index>` (all `_1` in this case since there are no parallel branches).

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"id"` (e.g. `"BR_1_10_1"`) | row index + `F_BUS`/`T_BUS` | ✓ derived |
| `"ports": {"bus1", "bus2"}` | `F_BUS` (col 1), `T_BUS` (col 2) | ✓ identical |
| `"params": {"R"}` | `BR_R` (col 3) | ✓ identical |
| `"params": {"X"}` | `BR_X` (col 4) | ✓ identical |
| `"params": {"B"}` | `BR_B` (col 5, total charging) | ✓ identical |
| `"params": {"G"}` | — (always 0.0) | ✓ |
| *(absent)* | `RATE_A/B/C` (cols 6–8) | ✗ not in JSON |
| *(absent)* | `TAP` (col 9), `SHIFT` (col 10) | ✗ not in JSON |
| *(absent)* | `BR_STATUS` (col 11) | ✗ not in JSON |

## Genrou generators (40 of 49 from .m)

**Total Genrou devices:** 40  
**Unique buses with generators:** 40 (one generator per bus — no multi-gen buses)

### Discrepancy vs. .m file

The `.m` file has **49 generators** in `mpc.gen`; `illinois.json` has **40 Genrou** devices. The discrepancy breaks down as follows:

- **11 offline** generators in `.m` (`GEN_STATUS=0`): buses 78, 79, 92, 161, 164, 165, 166, 168, 169, 196, 197
- **9 of those 11 are absent from JSON**: buses 78, 79, 92, 164, 165, 166, 168, 169, 196 (all `ng` fuel)
- **2 offline gens ARE included in JSON** despite `GEN_STATUS=0`: buses **161** (SPRINGFIELD 4 1, 138.6 MW) and **197** (GIBSON CITY 1 2, 67.5 MW)

This differs from the Hawaii case where **all** offline generators are excluded. Buses 161 and 197 appear to have been retained in the JSON dynamics model (possibly representing machines with dynamics that still need to be present in the simulation even when offline).

Net: 49 − 9 excluded = **40 Genrou** in JSON. ✓

All 40 Genrou devices have one generator per bus (no fanout needed for visualization).

### Field mapping: Genrou params vs. .m file

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"id"` (e.g. `"49_1"`) | `GEN_BUS` + rank within bus group in `mpc.gen` | ✓ derived; all Illinois gens are 1-per-bus, so all IDs end in `_1` |
| `"ports": {"bus"}` | `GEN_BUS` (col 1 of `mpc.gen`) | ✓ identical |
| `"params": {"p0"}` | `PG` (col 2) ÷ `mva` (machine MVA base) | ✓ per-unit on machine base |
| `"params": {"q0"}` | `QG` (col 3) ÷ `mva` | ✓ per-unit on machine base |
| `"params": {"mva"}` | `MBASE` (col 7, machine MVA base) | ✓ identical |
| `"params": {"H", "D", "Ra", "Xd", …}` | `.dyr` dynamics file (not in `.m`) | ✓ dynamics params from PowerWorld `.dyr` export |
| `"mon": ["delta", "omega"]` | — (GridKit-specific monitoring config) | n/a |
| *(absent)* | `QMAX`/`QMIN` (cols 4–5) | ✗ not in JSON |
| *(absent)* | `VG` (col 6, voltage setpoint) | ✗ not in JSON |
| *(absent)* | `PMAX`/`PMIN` (cols 9–10) | ✗ not in JSON |
| *(absent)* | `GEN_STATUS` (col 8) | ✗ not in JSON — most offline gens are simply absent; 2 exceptions (buses 161, 197) are included with their offline `p0`/`q0` = 0 |

### Excluded generators

9 offline generators excluded from JSON (all `ng` fuel):

| gen_row | bus | bus_name        | PMAX (MW) | JSON id (absent) |
|--------:|----:|:----------------|----------:|:-----------------|
|  16 |  78 | BRIMFIELD 3     |  18.0 | `78_1`  |
|  17 |  79 | BRIMFIELD 4     |  18.0 | `79_1`  |
|  20 |  92 | CLINTON 3 3     |   6.3 | `92_1`  |
|  38 | 164 | CHAMPAIGN 1 1   |  12.0 | `164_1` |
|  39 | 165 | CHAMPAIGN 1 2   |  26.0 | `165_1` |
|  40 | 166 | CHAMPAIGN 1 3   |   9.4 | `166_1` |
|  42 | 168 | CHAMPAIGN 1 5   |   9.4 | `168_1` |
|  43 | 169 | CHAMPAIGN 1 6   |   9.4 | `169_1` |
|  48 | 196 | GIBSON CITY 1 1 |  67.5 | `196_1` |

### Included buses with generators

All 40 Genrou buses (one generator per bus):

| # gens | Buses |
|-------:|:------|
| 1 each | 49, 50, 51, 52, 53, 65, 67, 68, 69, 70, 71, 72, 73, 76, 77, 90, 91, 94, 104, 105, 114, 115, 125, 126, 127, 135, 136, 147, 151, 152, 153, 154, 155, 161†, 167, 170, 182, 183, 189, 197† |

† Offline in `.m` (`GEN_STATUS=0`) but still present in JSON.

## BusFault devices (200 — parametric fault)

200 `BusFault` devices, one per bus (IDs `fault_1` through `fault_200`, ports `bus` 2–200 approximately). All have `"params": {"state0": false, "R": 0.0, "X": 0.001}`.

Unlike the **Hawaii** case — which has a single `BusFault` at bus 1 (ALOHA138) with `state0: true` (fault is active at $t=0$) — the Illinois case has faults at **every bus** but all with `state0: false` (disabled by default). The fault location is a **study variable** controlled externally per simulation run rather than being hard-coded in the model.

## Device model triplet: Genrou + Tgov1 + SexsPti

Every generator in the Illinois case is modeled as a triplet of three coupled devices:

- **Genrou** — synchronous machine (rotor electromechanics). Tracks rotor angle `delta`, rotor speed `omega`, and internal flux linkages. Inertia `H` is the epistemic uncertainty parameter for UQ campaigns.
- **Tgov1** — simplified steam/gas turbine governor. Regulates mechanical power `pmech` in response to frequency deviations.
- **SexsPti** — PTI SEXS simplified exciter/AVR (replaces the `Ieeet1` used in Hawaii). Regulates field voltage `efd` in response to terminal voltage deviations. Has 6 parameters: `Ta`, `Tb`, `Te`, `K`, `Efdmax`, `Efdmin`.

The exciter model difference vs Hawaii: Illinois uses **SexsPti** (simplified SEXS model with lead-lag filter and gain limiting) rather than **Ieeet1** (IEEE Type 1 exciter with more parameters and a transient gain reduction filter). Both regulate `efd` into the Genrou machine.

Loads are modeled as **LoadZIP** (164 devices), a ZIP load model mixing constant impedance (Z), constant current (I), and constant power (P) components — distinct from the simple `Load` (constant power) model used in Hawaii.

### Example JSON pieces (bus 49, id `49_1`)

```json
{
    "class": "Genrou",
    "ports": {"bus": 49, "speed": 0, "pmech": 1, "efd": 2},
    "id": "49_1",
    "params": {
        "p0": 0.04530000314116478,
        "q0": 0.0210999995470047,
        "H": 5.21145296,
        "D": 0.0,
        "Ra": 0.0,
        "Tdop": 6.81640005,
        "Tdopp": 0.0245,
        "Tqop": 1.39760005,
        "Tqopp": 0.0796,
        "Xd": 2.49710011,
        "Xdp": 0.43959999,
        "Xdpp": 0.31650001,
        "Xq": 2.36619997,
        "Xqp": 0.76020002,
        "Xqpp": 0.31650001,
        "Xl": 0.2384,
        "S10": 0.1971,
        "S12": 0.62940001,
        "mva": 5.44000006
    }
}
```

```json
{
    "class": "SexsPti",
    "ports": {"bus": 49, "efd": 2},
    "id": "49_1_sexs_pti",
    "params": {
        "Ta": 0.80000001,
        "Tb": 8.0,
        "Te": 0.05,
        "K": 130.0,
        "Efdmax": 5.0,
        "Efdmin": -4.0
    }
}
```
