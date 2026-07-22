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

This differs from the Hawaii case where **all** offline generators are excluded. Buses 161 and 197 have `GEN_STATUS=0` in this particular `.m` solution file, but their Genrou devices in the JSON have non-zero p0/q0. This suggests the JSON was built from a different snapshot of the `.m` case (likely the base case) where those generators were online, and the JSON was not updated to reflect the offline state in subsequent solutions.

Net: 49 − 9 excluded = **40 Genrou** in JSON. ✓

All 40 Genrou devices have one generator per bus (no fanout needed for visualization).

### Field mapping: Genrou params vs. .m file

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"id"` (e.g. `"49_1"`) | `GEN_BUS` + rank within bus group in `mpc.gen` | ✓ derived; all Illinois gens are 1-per-bus, so all IDs end in `_1` |
| `"ports": {"bus"}` | `GEN_BUS` (col 1 of `mpc.gen`) | ✓ identical |
| `"params": {"p0"}` | `PG` (col 2) ÷ `baseMVA` (100 MVA) — conversion assumed; see table below | see table |
| `"params": {"q0"}` | `QG` (col 3) ÷ `baseMVA` — conversion assumed; see table below | see table |
| `"params": {"mva"}` | `MBASE` (col 7, machine MVA base) | ✓ identical |
| `"params": {"H", "D", "Ra", "Xd", …}` | `.dyr` dynamics file (not in `.m`) | ✓ dynamics params from PowerWorld `.dyr` export |
| `"mon": ["delta", "omega"]` | — (GridKit-specific monitoring config) | n/a |
| *(absent)* | `QMAX`/`QMIN` (cols 4–5) | ✗ not in JSON |
| *(absent)* | `VG` (col 6, voltage setpoint) | ✗ not in JSON |
| *(absent)* | `PMAX`/`PMIN` (cols 9–10) | ✗ not in JSON |
| *(absent)* | `GEN_STATUS` (col 8) | ✗ not in JSON — most offline gens are simply absent; 2 exceptions (buses 161, 197) are included with non-zero p0/q0. These values reflect the operating point used to build the JSON, not the TAMU base case `.m`. |

See [Patching a PF solution — Generator dispatch](#generator-dispatch-pgqg--genrou-p0q0) for the full per-generator comparison table.

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

## Patching a PF solution (.m) into the JSON case

To transfer a power flow solution from an hourly `.m` file into `illinois.json`, three sections of the JSON must be updated. Everything else (dynamics params, topology, fault devices) is static.

### What to patch

| JSON section | Field(s) patched | Source in `.m` file | Formula |
|---|---|---|---|
| `buses[i].init` | `Vr`, `Vi` | `mpc.bus` col 8 `VM` (pu), col 9 `VA` (degrees) | `Vr = VM × cos(VA × π/180)`, `Vi = VM × sin(VA × π/180)` |
| `Genrou.params` (online gens) | `p0`, `q0` | `mpc.gen` col 2 `PG` (MW), col 3 `QG` (MVAr) | `p0 = PG / 100`, `q0 = QG / 100` |
| `LoadZIP.params` | `Pnom`, `Qnom`, `Vnom` | `mpc.bus` col 3 `PD` (MW), col 4 `QD` (MVAr), col 8 `VM` (pu) | `Pnom = PD / 100`, `Qnom = QD / 100`, `Vnom = VM` |

Notes:
- `va_base = 1e8` VA = 100 MVA (from `header.va_base`), so dividing MW by 100 gives pu on system base.
- `Vnom` in LoadZIP is the per-unit bus voltage magnitude; it serves as the ZIP model reference voltage.
- `Tgov1`, `SexsPti`, and `Branch` params are pure dynamics/topology — never change between hours.
- `BusFault.state0` is scenario-controlled (which bus to fault), not derived from the PF solution.

In all three cases the current JSON values do **not** match the TAMU base case `.m` file — they were derived from a different (likely PowerWorld) power flow solution. The formulas above are still the correct conversion; they will produce values consistent with any `.m` solution file that has been run through MATPOWER AC power flow.

---

### Bus voltages: VM/VA → Vr/Vi

Formula: `Vr = VM × cos(VA × π/180)`, `Vi = VM × sin(VA × π/180)`

Sample comparison of the first 10 buses between the TAMU base case `.m` and the JSON init values. `Vr_calc` / `Vi_calc` are what the formula produces from the `.m` file; `Vr (json)` / `Vi (json)` are the current JSON values (different PF solution).

| bus | VM (.m) | VA° (.m) | Vr_calc | Vi_calc | Vr (json) | Vi (json) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.019152 | -7.0852 | 1.011370 | -0.125707 | 0.230636 | -0.973007 |
| 2 | 1.019035 | -7.0980 | 1.011225 | -0.125919 | 0.230264 | -0.972914 |
| 3 | 1.030054 | -10.0291 | 1.014314 | -0.179382 | 0.212164 | -0.995483 |
| 4 | 1.030029 | -10.0320 | 1.014280 | -0.179429 | 0.212075 | -0.995461 |
| 5 | 1.037273 | -3.5463 | 1.035287 | -0.064160 | 0.144770 | -1.016041 |
| 6 | 1.037147 | -3.5596 | 1.035146 | -0.064392 | 0.144417 | -1.015916 |
| 7 | 1.013098 | -7.7932 | 1.003741 | -0.137374 | 0.208815 | -0.968258 |
| 8 | 1.012736 | -7.8331 | 1.003287 | -0.138023 | 0.207668 | -0.967950 |
| 9 | 1.015774 | -7.4718 | 1.007149 | -0.132089 | 0.220768 | -0.970143 |
| 10 | 1.015051 | -7.5619 | 1.006224 | -0.133578 | 0.218199 | -0.969615 |

**Summary**: All 200 buses are present in both `.m` and JSON. The JSON |V| magnitudes range 0.986–1.040 pu (mean 1.020); the TAMU base case VM range is 1.010–1.056 pu (mean 1.034). Angles differ substantially — the JSON was built from a completely different load flow. After patching, Vr/Vi will reflect the hourly PF solution.

---

### Generator dispatch: PG/QG → Genrou p0/q0

Formula: `p0 = PG / 100`, `q0 = QG / 100` (pu on 100 MVA system base)

All 49 generators from the TAMU base case `.m`, compared to the current JSON values. `p0 match` / `q0 match` = within 1e-3. Rows with *(absent)* are excluded from the JSON entirely and are never patched.

| gen_row | bus | bus_name | .m status | JSON id | p0 (json) | q0 (json) | PG/100 | QG/100 | p0 match | q0 match |
|---:|---:|:---|:---:|:---:|---:|---:|---:|---:|:---:|:---:|
| 1 | 49 | RANTOUL 2 1 | on | `49_1` | 0.0453 | 0.0211 | 0.0136 | 0.0088 | ✗ | ✗ |
| 2 | 50 | RANTOUL 2 2 | on | `50_1` | 0.0453 | 0.0211 | 0.0136 | 0.0118 | ✗ | ✗ |
| 3 | 51 | RANTOUL 2 3 | on | `51_1` | 0.0453 | 0.0211 | 0.0136 | 0.0077 | ✗ | ✗ |
| 4 | 52 | RANTOUL 2 4 | on | `52_1` | 0.0453 | 0.0211 | 0.0136 | 0.0125 | ✗ | ✗ |
| 5 | 53 | RANTOUL 2 5 | on | `53_1` | 0.0907 | 0.0248 | 0.0272 | 0.0079 | ✗ | ✗ |
| 6 | 65 | PAXTON 1 1 | on | `65_1` | 0.0400 | 0.1101 | 0.8650 | 0.0084 | ✗ | ✗ |
| 7 | 67 | MOUNT ZION 1 | on | `67_1` | 0.0470 | 0.0168 | 0.0141 | -0.0057 | ✗ | ✗ |
| 8 | 68 | MOUNT ZION 2 | on | `68_1` | 0.2792 | 0.0151 | 0.0838 | -0.0055 | ✗ | ✗ |
| 9 | 69 | MOUNT ZION 3 | on | `69_1` | 0.2792 | 0.0263 | 0.0838 | -0.0153 | ✗ | ✗ |
| 10 | 70 | MOUNT ZION 4 | on | `70_1` | 0.2792 | 0.0146 | 0.0838 | -0.0078 | ✗ | ✗ |
| 11 | 71 | MOUNT ZION 5 | on | `71_1` | 0.2792 | 0.0326 | 0.0838 | -0.0170 | ✗ | ✗ |
| 12 | 72 | MOUNT ZION 6 | on | `72_1` | 0.2792 | 0.0351 | 0.0838 | -0.0186 | ✗ | ✗ |
| 13 | 73 | MOUNT ZION 7 | on | `73_1` | 0.2792 | 0.0154 | 0.0838 | -0.0111 | ✗ | ✗ |
| 14 | 76 | BRIMFIELD 1 | on | `76_1` | 0.0120 | 0.0204 | 0.0120 | 0.0141 | ✓ | ✗ |
| 15 | 77 | BRIMFIELD 2 | on | `77_1` | 0.0072 | 0.0116 | 0.0072 | 0.0050 | ✓ | ✗ |
| 16 | 78 | BRIMFIELD 3 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 17 | 79 | BRIMFIELD 4 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 18 | 90 | CLINTON 3 1 | on | `90_1` | 0.0096 | 0.0108 | 0.0096 | 0.0013 | ✓ | ✗ |
| 19 | 91 | CLINTON 3 2 | on | `91_1` | 0.0150 | 0.0255 | 0.0150 | 0.0059 | ✓ | ✗ |
| 20 | 92 | CLINTON 3 3 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 21 | 94 | TUSCOLA 2 1 | on | `94_1` | 0.1800 | 0.0163 | 0.0540 | -0.0220 | ✗ | ✗ |
| 22 | 104 | ELLSWORTH 1 2 | on | `104_1` | 0.1790 | 0.1549 | 0.6760 | 0.0310 | ✗ | ✗ |
| 23 | 105 | ELLSWORTH 1 3 | on | `105_1` | 0.3360 | 0.2308 | 1.5480 | 0.0829 | ✗ | ✗ |
| 24 | 114 | NORMAL 2 2 | on | `114_1` | 0.0040 | 0.0036 | 0.0140 | 0.0036 | ✗ | ✓ |
| 25 | 115 | NORMAL 2 3 | on | `115_1` | 0.4770 | 0.3099 | 1.3350 | 0.1260 | ✗ | ✗ |
| 26 | 125 | BARTONVILLE 2 | on | `125_1` | 1.5000 | 0.2195 | 0.3902 | 0.0724 | ✗ | ✗ |
| 27 | 126 | BARTONVILLE 3 | on | `126_1` | 1.5000 | 0.3111 | 0.3902 | 0.1255 | ✗ | ✗ |
| 28 | 127 | BARTONVILLE 4 | on | `127_1` | 1.5000 | 0.5593 | 0.3902 | 0.2599 | ✗ | ✗ |
| 29 | 135 | PEKIN 1 2 | on | `135_1` | 1.9642 | 0.5162 | 1.3392 | 0.2146 | ✗ | ✗ |
| 30 | 136 | PEKIN 1 3 | on | `136_1` | 2.2255 | 0.6153 | 1.3392 | 0.2506 | ✗ | ✗ |
| 31 | 147 | HOPEDALE 2 1 | on | `147_1` | 0.1340 | 0.2141 | 0.9240 | 0.0879 | ✗ | ✗ |
| 32 | 151 | SPRINGFIELD 5 2 | on | `151_1` | 0.0540 | 0.0105 | 0.0162 | 0.0002 | ✗ | ✗ |
| 33 | 152 | SPRINGFIELD 5 3 | on | `152_1` | 0.7722 | 0.0770 | 0.2317 | 0.0027 | ✗ | ✗ |
| 34 | 153 | SPRINGFIELD 5 4 | on | `153_1` | 0.7722 | 0.0676 | 0.2317 | 0.0002 | ✗ | ✗ |
| 35 | 154 | SPRINGFIELD 5 5 | on | `154_1` | 0.7722 | 0.0804 | 0.2317 | 0.0014 | ✗ | ✗ |
| 36 | 155 | SPRINGFIELD 5 6 | on | `155_1` | 0.7722 | 0.0812 | 0.2317 | 0.0013 | ✗ | ✗ |
| 37 | 161 | SPRINGFIELD 4 1 | **off** | `161_1`† | 0.4158 | 0.1800 | 0.0000 | 0.0000 | ✗ | ✗ |
| 38 | 164 | CHAMPAIGN 1 1 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 39 | 165 | CHAMPAIGN 1 2 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 40 | 166 | CHAMPAIGN 1 3 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 41 | 167 | CHAMPAIGN 1 4 | on | `167_1` | 0.0282 | 0.0478 | 0.0282 | -0.0104 | ✓ | ✗ |
| 42 | 168 | CHAMPAIGN 1 5 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 43 | 169 | CHAMPAIGN 1 6 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 44 | 170 | CHAMPAIGN 1 7 | on | `170_1` | 0.0282 | 0.0478 | 0.0282 | -0.0099 | ✓ | ✗ |
| 45 | 182 | SPRINGFIELD 2 1 | on | `182_1` | 0.1750 | 0.0498 | 0.0525 | 0.0215 | ✗ | ✗ |
| 46 | 183 | SPRINGFIELD 2 2 | on | `183_1` | 0.2660 | 0.1240 | 0.0798 | 0.0573 | ✗ | ✗ |
| 47 | 189 | CLINTON 1 2 | on | `189_1` | 5.6914 | 0.4956 | 3.8437 | -0.2415 | ✗ | ✗ |
| 48 | 196 | GIBSON CITY 1 1 | **off** | *(absent)* | — | — | 0.0000 | 0.0000 | — | — |
| 49 | 197 | GIBSON CITY 1 2 | **off** | `197_1`† | 0.2025 | 0.0902 | 0.0000 | 0.0000 | ✗ | ✗ |

† Present in JSON despite `GEN_STATUS=0` in the TAMU base case `.m`; see [Offline generator handling](#offline-generator-handling).

**Summary**: p0/q0 in the JSON almost never match the TAMU base case `.m`. Only 6 of 40 online generators have p0 within 1e-3 (rows 14, 15, 18, 19, 41, 44 — all small ng units). The JSON p0/q0 were derived from a different PF solution (likely PowerWorld). After patching with `PG / 100` from each hourly `.m` solution, values will reflect the actual dispatch.

---

### Load demand: PD/QD/VM → LoadZIP Pnom/Qnom/Vnom

Formula: `Pnom = PD / 100`, `Qnom = QD / 100`, `Vnom = VM`

**Coverage**: 108 buses have `PD > 0` in the TAMU base case `.m`; all 108 have a corresponding LoadZIP in the JSON. One additional bus (bus 15) has a LoadZIP in the JSON but `PD = 0` in the base case `.m` — see extra-bus note below. Total LoadZIP devices: 164 across 109 buses; 38 buses have 2 LoadZIP devices each (see multi-ZIP examples below).

Sample comparison, base case `.m` vs current JSON. Rows marked `(×2)` have 2 LoadZIP devices; individual device values are shown indented:

| bus | PD (MW) | QD (MVAr) | VM (.m) | PD/100 | QD/100 | #ZIP | sum_Pnom (json) | sum_Qnom (json) | Vnom (json) | Pnom match | Qnom match |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 2 | 7.39 | 2.10 | 1.0190 | 0.0739 | 0.0210 | 1 | 0.1082 | 0.0308 | 0.9998 | ✗ | ✗ |
| 4 | 1.70 | 0.48 | 1.0300 | 0.0170 | 0.0048 | 1 | 0.0268 | 0.0076 | 1.0178 | ✗ | ✗ |
| 6 | 7.95 | 2.26 | 1.0371 | 0.0795 | 0.0226 | 1 | 0.1082 | 0.0308 | 1.0261 | ✗ | ✗ |
| 8 | 23.74 | 6.77 | 1.0127 | 0.2374 | 0.0677 | 1 | 0.3477 | 0.0991 | 0.9900 | ✗ | ✗ |
| 10 | 42.68 | 12.16 | 1.0151 | 0.4268 | 0.1216 | 1 | 0.6253 | 0.1782 | 0.9939 | ✗ | ✗ |
| **39** (×2) | **7.41** | **2.11** | **1.0246** | **0.0741** | **0.0211** | **2** | **0.1172** | **0.0334** | **1.0088** | ✗ | ✗ |
| ↳ ZIP[0] | | | | | | | 0.0162 | 0.0046 | 1.0088 | | |
| ↳ ZIP[1] | | | | | | | 0.1010 | 0.0288 | 1.0088 | | |
| **42** (×2) | **10.62** | **3.03** | **1.0301** | **0.1062** | **0.0303** | **2** | **0.1681** | **0.0479** | **1.0191** | ✗ | ✗ |
| ↳ ZIP[0] | | | | | | | 0.0169 | 0.0048 | 1.0191 | | |
| ↳ ZIP[1] | | | | | | | 0.1512 | 0.0431 | 1.0191 | | |
| **44** (×2) | **4.15** | **1.18** | **1.0280** | **0.0415** | **0.0118** | **2** | **0.0657** | **0.0187** | **1.0150** | ✗ | ✗ |
| ↳ ZIP[0] | | | | | | | 0.0150 | 0.0043 | 1.0150 | | |
| ↳ ZIP[1] | | | | | | | 0.0506 | 0.0144 | 1.0150 | | |
| 16 | 44.76 | 12.76 | 1.0485 | 0.4476 | 0.1276 | 1 | 0.7086 | 0.2019 | 1.0371 | ✗ | ✗ |
| 18 | 0.71 | 0.20 | 1.0248 | 0.0071 | 0.0020 | 1 | 0.0113 | 0.0032 | 1.0059 | ✗ | ✗ |
| **15** *(extra)* | **0** | **0** | **1.0491** | **0** | **0** | **1** | **0.0000** | **-0.3233** | **1.0381** | n/a | n/a |

**Bus 15 (extra LoadZIP)**: `PD=QD=0` in the TAMU base case `.m`, but the JSON has one LoadZIP with `Pnom=0, Qnom=-0.3233`. Negative Qnom indicates a **capacitor bank** (reactive compensation). This device is not a real load and should not be patched with `PD/QD` from the `.m` file; its Qnom should remain as-is (or be determined by a separate reactive dispatch).

**Multi-ZIP buses**: When a bus has 2 LoadZIP devices, each represents a separate load at that bus (e.g., different customers or load classes). The sum of all `Pnom`/`Qnom` values across devices at the same bus is what must match `PD/100` / `QD/100` after patching. All devices share the same `Vnom` (the bus voltage from the PF solution). In `m_to_case.py`, the total `PD` from the `.m` file must be distributed across the N devices at each bus — the simplest approach is to scale each device's existing `Pnom` by the ratio `(PD/100) / sum_Pnom_json`.

**Summary**: JSON `Pnom` values are systematically ~1.49× higher than `PD/100` from the TAMU base case (ratio range 1.36–1.60 across all buses). This is a higher-load operating point used to build the JSON, not a fixed scale factor. After patching, `Pnom`/`Qnom` will reflect the actual hourly load from each `.m` solution.

---

### Branches: nothing to patch

`Branch` devices in the JSON store only static network data: `R`, `X`, `G`, `B` (pi-model impedance and shunt admittance) plus `ports.bus1` / `ports.bus2`. These are derived from the `.m` branch table and do not change between hours. Branch flows (PF, QF, PT, QT) are outputs of the transient simulation, not inputs, so there is nothing to write from the PF solution.

`RATE_A` from the `.m` branch table is used only in `m_viz_utils.py` for the loading-percentage color scale in the geo map; it does not appear in the JSON and is not needed by the transient solver.

**Transformer tap ratios**: `case_ACTIVSg200.m` has 245 branches, all with `ratio=0` (179 lines, MATPOWER convention for "treat as 1") or `ratio=1.0` (66 transformer entries). There are zero non-unity off-nominal taps in this case, so no special tap-ratio handling is needed in `m_to_case.py` or in the GridKit PF solver.

---

### Offline generator handling

The JSON has no mechanism to disable a Genrou at runtime. The rule is:

| Generator state in `.m` | Generator in JSON | Action in `m_to_case.py` |
|---|---|---|
| `GEN_STATUS=1` (online) | present | patch `p0`, `q0` |
| `GEN_STATUS=0` (offline) | absent (9 gens: buses 78, 79, 92, 164–166, 168, 169, 196) | nothing (already absent) |
| `GEN_STATUS=0` (offline) | present (buses 161, 197) | open question — remove device or keep with base-case p0/q0? |

The two anomalous gens (buses 161, 197) are in the JSON with non-zero p0/q0 from the original build operating point. Whether to keep or remove them when they are offline in the hourly `.m` solution needs to be decided in consultation with the GridKit dev team.

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

---

## UQ parameter selection

### Fault location

Fault is applied at **bus 2**, using the first `BusFault` device in `illinois.json` (`element_id=0`). The base `illinois.solver.json` has tmax=20.0 s with the fault at t=10.0/10.1 s. For v1 we override to tmax=10.0 s and fault at t=1.0/1.1 s (matching the Hawaii cadence) to avoid running 10 s of pre-fault steady state.

### Generator selection

Generators were selected at increasing hop distance from the fault bus (BFS on the branch adjacency graph). One generator per hop level gives a natural set spanning close-in to far-out.

| id | bus | bus_name | hop distance | fuel | PMAX (MW) | H (nominal) |
|----|----:|:---------|:------------:|:----:|----------:|------------:|
| `126_1` | 126 | BARTONVILLE 3  | 4 | coal    | 130.05 | 7.29971 |
| `135_1` | 135 | PEKIN 1 2      | 5 | coal    | 446.40 | 3.23290 |
| `115_1` | 115 | NORMAL 2 3     | 7 | wind    | 133.50 | 2.72678 |
| `189_1` | 189 | CLINTON 1 2    | 9 | nuclear | 569.15 | 3.40165 |

All four are Genrou machines with varying H values. Since H is the UQ parameter here, the spread of nominal H values across the selected generators (2.7–7.3 s) is the primary source of diversity in the ensemble.

### v1 (1K samples) — epistemic UQ, 4 generators

- dist: **Gaussian**, std = 12% of nominal (epistemic uncertainty in H)
- method: LHS + normal ppf (scipy), seed=42
- serialize: per-run (`runs/run_NNN.parquet`)
- fault: bus 2, t=1.0 s on / t=1.1 s cleared, tmax=10.0 s (solver override; base has tmax=20 s, fault at t=10/10.1 s)
- data path: `/kfs2/projects/scidac/scidac-data/gridkit-runs/illinois-v1/`
- metadata: `illinois-v1/meta.yml`, samples: `illinois-v1/samples.csv`

| id | bus | H nominal | mean | std (12%) |
|----|----:|----------:|-----:|----------:|
| `126_1` | 126 | 7.29971 | 7.29971 | 0.87597 |
| `135_1` | 135 | 3.23290 | 3.23290 | 0.38795 |
| `115_1` | 115 | 2.72678 | 2.72678 | 0.32721 |
| `189_1` | 189 | 3.40165 | 3.40165 | 0.40820 |

### v2 (4K samples) — epistemic UQ, 4 generators

Same as v1 with N increased to 4000 for better coverage of the 4D input space.

- dist: **Gaussian**, std = 12% of nominal (epistemic uncertainty in H)
- method: LHS + normal ppf (scipy), seed=42
- serialize: per-run (`runs/run_NNN.parquet`)
- fault: bus 2, t=1.0 s on / t=1.1 s cleared, tmax=10.0 s (solver override)
- data path: `/kfs2/projects/scidac/scidac-data/gridkit-runs/illinois-v2/`
- metadata: `illinois-v2/meta.yml`, samples: `illinois-v2/samples.csv`

| id | bus | H nominal | mean | std (12%) |
|----|----:|----------:|-----:|----------:|
| `126_1` | 126 | 7.29971 | 7.29971 | 0.87597 |
| `135_1` | 135 | 3.23290 | 3.23290 | 0.38795 |
| `115_1` | 115 | 2.72678 | 2.72678 | 0.32721 |
| `189_1` | 189 | 3.40165 | 3.40165 | 0.40820 |
