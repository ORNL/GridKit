# Hawaii case — MATPOWER .m file

## Geo plot

Interactive Plotly map — branch loading (viridis), buses sized by PD, generators, fault marker at ALOHA138.  
Generated from `m_viz.ipynb` with `show_loading=True`, `fault_bus="ALOHA138"`.

[Open map in new tab (larger view) ↗](Hawaii40_geo_loading_faultALOHA138.html){target="_blank"}

<iframe src="Hawaii40_geo_loading_faultALOHA138.html" width="100%" height="650px" style="border:none;"></iframe>

## Bus table (37 buses)

BUS_TYPE: 1 = PQ load bus, 2 = PV generator bus, 3 = reference/slack bus.  
`BASE_KV` is the voltage level in kV — the numbers in bus names (e.g. `138`, `69`) refer to this, **not** the bus index.

| BUS_I | bus_name      | BUS_TYPE | BASE_KV | PD (MW) | QD (MVAr) |
|------:|:--------------|:--------:|--------:|--------:|----------:|
|  1    | ALOHA138      | 1        | 138     |   0.00  |    0.0    |
|  2    | ALOHA69       | 2 (PV)   |  69     |  60.70  |    0.0    |
|  3    | FLOWER69      | 1        |  69     |  59.39  |    0.0    |
|  4    | WAVE69        | 1        |  69     |  22.47  |    0.0    |
|  5    | HONOLULU138   | 1        | 138     |   0.00  |    0.0    |
|  6    | HONOLULU69    | 1        |  69     |  27.46  |    0.0    |
|  7    | SURF69        | 1        |  69     |  37.01  |    0.0    |
|  8    | KANEOHE69     | 1        |  69     |  65.10  |    0.0    |
|  9    | TURTLE138     | 1        | 138     |   0.00  |    0.0    |
| 10    | TURTLE69      | 1        |  69     |  37.76  |    0.0    |
| 11    | MAHALO69      | 1        |  69     |  65.55  |    0.0    |
| 12    | LYCHEE69      | 1        |  69     |  54.01  |    0.0    |
| 13    | COCONUT69     | 1        |  69     |  59.24  |    0.0    |
| 14    | KAILUA138     | 1        | 138     |   0.00  |    0.0    |
| 15    | KAILUA69      | 1        |  69     |  60.90  |    0.0    |
| 16    | PALM69        | 1        |  69     |  23.83  |    0.0    |
| 17    | WAIMANALO69   | 1        |  69     |  12.04  |    0.0    |
| 18    | VOLCANO69     | 1        |  69     |  36.32  |    0.0    |
| 19    | PEARL CITY69  | 1        |  69     |  48.60  |    0.0    |
| 20    | MILILANI69    | 1        |  69     |  64.95  |    0.0    |
| 21    | AIEA69        | 1        |  69     |  48.34  |    0.0    |
| 22    | WAIPAHU138    | 1        | 138     |   0.00  |    0.0    |
| 23    | WAIPAHU69     | 3 (ref)  |  69     |  86.75  |    0.0    |
| 24    | KAPOLEI69     | 1        |  69     |  46.58  |    0.0    |
| 25    | EWA BEACH138  | 1        | 138     |   0.00  |    0.0    |
| 26    | EWA BEACH69   | 2 (PV)   |  69     |  75.28  |    0.0    |
| 27    | KAHUKU69      | 2 (PV)   |  69     |   3.95  |    0.0    |
| 28    | HALEIWA69     | 2 (PV)   |  69     |   8.82  |    0.0    |
| 29    | LAIE69        | 1        |  69     |   7.70  |    0.0    |
| 30    | WAHIAWA69     | 1        |  69     |  49.03  |    0.0    |
| 31    | WAIALUA69     | 1        |  69     |   9.62  |    0.0    |
| 32    | HAUULA69      | 1        |  69     |   6.67  |    0.0    |
| 33    | WAIANAE69     | 2 (PV)   |  69     |  58.22  |    0.0    |
| 34    | SCHOFIELD69   | 2 (PV)   |  69     |   0.00  |    0.0    |
| 35    | KALAELOA138   | 2 (PV)   | 138     |   0.00  |    0.0    |
| 36    | COGEN69       | 2 (PV)   |  69     |   0.00  |    0.0    |
| 37    | KAHE138       | 2 (PV)   | 138     |   0.00  |    0.0    |

## Gen table (45 generators)

`gen_row` is the 1-based row index in the MATPOWER `mpc.gen` matrix (same as `GEN_BUS` device index in GridKit). `PG`/`QG` are dispatch set-points; `PMAX`/`PMIN` are capacity limits — all in MW or MVAr.

| gen_row | GEN_BUS | bus_name     | PG (MW) | QG (MVAr) | PMAX  | PMIN |
|--------:|--------:|:-------------|--------:|----------:|------:|-----:|
|  1–4    |  2      | ALOHA69      |   2.50  |    0.80   |   2.5 |  1.0 |
|  5      | 23      | WAIPAHU69    |  69.27  |    0.04   |  77.8 | 20.0 |
|  6      | 23      | WAIPAHU69    |  39.13  |    0.04   |  45.9 |  0.0 |
|  7      | 23      | WAIPAHU69    |  69.27  |    0.04   |  77.8 | 20.0 |
|  8      | 23      | WAIPAHU69    |  47.64  |    0.04   |  51.9 | 23.0 |
|  9      | 23      | WAIPAHU69    |  47.11  |    0.04   |  51.8 | 20.0 |
| 10      | 23      | WAIPAHU69    |  44.53  |    0.04   |  51.2 |  6.0 |
| 11      | 23      | WAIPAHU69    |  44.53  |    0.04   |  51.2 |  6.0 |
| 12      | 23      | WAIPAHU69    |  43.61  |    0.04   |  47.7 | 20.0 |
| 13      | 23      | WAIPAHU69    |  43.19  |    0.04   |  47.2 | 20.0 |
| 14      | 23      | WAIPAHU69    |  12.53  |    0.04   |  14.7 |  0.0 |
| 15      | 26      | EWA BEACH69  |  20.00  |    6.00   |  20.0 |  0.0 |
| 16      | 26      | EWA BEACH69  |  10.20  |    3.10   |  10.2 |  0.0 |
| 17      | 27      | KAHUKU69     |  30.00  |   -5.40   |  30.0 |  0.0 |
| 18      | 27      | KAHUKU69     |  27.60  |   -5.00   |  27.6 |  1.0 |
| 19      | 28      | HALEIWA69    |  46.30  |   -5.49   |  49.0 |  0.0 |
| 20      | 28      | HALEIWA69    |  69.00  |   -5.49   |  69.0 |  0.0 |
| 21      | 33      | WAIANAE69    |  27.60  |    8.30   |  27.6 |  0.0 |
| 22–27   | 34      | SCHOFIELD69  |  ~6–7   |   -0.98   |   8.4 |  2.5 |
| 28      | 35      | KALAELOA138  |  58.00  |    0.05   |  58.0 |  5.0 |
| 29      | 35      | KALAELOA138  |  28.00  |    0.05   |  28.0 |  5.0 |
| 30      | 35      | KALAELOA138  |   0.00  |    0.00   | 180.0 | 63.0 |
| 31      | 35      | KALAELOA138  |  20.00  |    0.05   |  20.0 |  5.0 |
| 32      | 35      | KALAELOA138  |  67.00  |    0.05   |  85.0 | 40.0 |
| 33      | 35      | KALAELOA138  |   0.00  |    0.00   |  85.0 | 40.0 |
| 34      | 35      | KALAELOA138  |  50.00  |    0.05   |  50.0 | 10.0 |
| 35      | 35      | KALAELOA138  |  56.20  |    0.05   | 113.0 | 42.0 |
| 36–39   | 36      | COGEN69      |  ~2–2.5 |    1.00   |   3.2 |  2.0 |
| 40      | 37      | KAHE138      |   0.00  |    0.00   | 128.7 | 45.0 |
| 41      | 37      | KAHE138      |   0.00  |    0.00   | 128.1 | 55.0 |
| 42      | 37      | KAHE138      |  58.32  |    5.63   |  87.2 | 20.0 |
| 43      | 37      | KAHE138      |   0.00  |    0.00   |  82.1 | 20.0 |
| 44      | 37      | KAHE138      |  68.48  |    5.63   |  78.1 | 30.0 |
| 45      | 37      | KAHE138      |   0.00  |    0.00   |  77.9 | 20.0 |

## Branch table (89 branches)

`TAP` = 0.0 means a transmission line (no transformer tap); `TAP` = 1.0 means an ideal transformer (ratio 1:1). Resistance and reactance are in per-unit on the system MVA base. `RATE_A` is the thermal rating in MVA.

| branch_row | F_BUS | T_BUS | F name       | T name       | BR_R   | BR_X    | RATE_A | TAP |
|-----------:|------:|------:|:-------------|:-------------|-------:|--------:|-------:|:---:|
|  1–3   |  1 | 2  | ALOHA138      | ALOHA69       | 0.00484 | 0.12377 | 123.2 | 1.0 (xfmr) |
|  4–5   |  1 | 5  | ALOHA138      | HONOLULU138   | 0.00340 | 0.02091 | 232.8 | —  |
|  6     |  1 | 22 | ALOHA138      | WAIPAHU138    | 0.00374 | 0.02417 | 177.8 | —  |
|  7–8   |  1 | 25 | ALOHA138      | EWA BEACH138  | 0.00289 | 0.01580 | 217.0 | —  |
|  9–10  |  1 | 35 | ALOHA138      | KALAELOA138   | 0.01080 | 0.04192 | 281.6 | —  |
| 11     |  2 |  4 | ALOHA69       | WAVE69        | 0.04296 | 0.07936 |  66.8 | —  |
| 12–13  |  2 |  6 | ALOHA69       | HONOLULU69    | 0.02689 | 0.05860 |  82.0 | —  |
| 14     |  2 |  7 | ALOHA69       | SURF69        | 0.02754 | 0.07569 |  77.8 | —  |
| 15     |  2 | 21 | ALOHA69       | AIEA69        | 0.03646 | 0.07734 |  94.8 | —  |
| 16–18  |  2 | 23 | ALOHA69       | WAIPAHU69     | 0.03262 | 0.08103 | 101.3 | —  |
| 19–20  |  2 | 26 | ALOHA69       | EWA BEACH69   | 0.03049 | 0.06612 |  61.7 | —  |
| 21–22  |  3 |  6 | FLOWER69      | HONOLULU69    | 0.01478 | 0.03669 |  85.0 | —  |
| 23     |  3 | 16 | FLOWER69      | PALM69        | 0.02383 | 0.06359 |  82.8 | —  |
| 24–26  |  3 | 19 | FLOWER69      | PEARL CITY69  | 0.01382 | 0.03123 |  91.6 | —  |
| 27     |  3 | 21 | FLOWER69      | AIEA69        | 0.02423 | 0.04751 |  49.0 | —  |
| 28     |  4 | 13 | WAVE69        | COCONUT69     | 0.01214 | 0.03849 |  53.2 | —  |
| 29     |  5 |  6 | HONOLULU138   | HONOLULU69    | 0.00479 | 0.13956 | 198.0 | 1.0 (xfmr) |
| 30     |  5 |  9 | HONOLULU138   | TURTLE138     | 0.00229 | 0.01206 | 170.7 | —  |
| 31     |  5 | 14 | HONOLULU138   | KAILUA138     | 0.00574 | 0.03882 | 162.5 | —  |
| 32–34  |  7 |  6 | SURF69        | HONOLULU69    | 0.00400 | 0.01000 |  40.2 | —  |
| 35     |  6 | 11 | HONOLULU69    | MAHALO69      | 0.01980 | 0.04240 |  89.0 | —  |
| 36–37  |  6 | 12 | HONOLULU69    | LYCHEE69      | 0.01364 | 0.03319 |  58.8 | —  |
| 38     |  7 | 10 | SURF69        | TURTLE69      | 0.00400 | 0.01000 |  50.0 | —  |
| 39–41  |  7 | 13 | SURF69        | COCONUT69     | 0.00941 | 0.02771 |  40.2 | —  |
| 42     |  8 | 15 | KANEOHE69     | KAILUA69      | 0.03139 | 0.07585 |  84.9 | —  |
| 43–44  |  8 | 19 | KANEOHE69     | PEARL CITY69  | 0.03699 | 0.10246 |  90.0 | —  |
| 45     |  9 | 10 | TURTLE138     | TURTLE69      | 0.00468 | 0.12565 |  56.2 | 1.0 (xfmr) |
| 46     |  9 | 14 | TURTLE138     | KAILUA138     | 0.00933 | 0.03785 | 169.3 | —  |
| 47     | 11 | 13 | MAHALO69      | COCONUT69     | 0.01844 | 0.04777 |  45.3 | —  |
| 48     | 11 | 15 | MAHALO69      | KAILUA69      | 0.02328 | 0.06922 |  48.6 | —  |
| 49     | 12 | 15 | LYCHEE69      | KAILUA69      | 0.01511 | 0.04146 |  77.0 | —  |
| 50     | 12 | 17 | LYCHEE69      | WAIMANALO69   | 0.02278 | 0.06534 |  49.9 | —  |
| 51     | 13 | 16 | COCONUT69     | PALM69        | 0.01309 | 0.03255 |  52.0 | —  |
| 52–53  | 14 | 15 | KAILUA138     | KAILUA69      | 0.00325 | 0.13152 |  85.9 | 1.0 (xfmr) |
| 54     | 16 | 18 | PALM69        | VOLCANO69     | 0.01553 | 0.03923 |  51.7 | —  |
| 55     | 17 | 18 | WAIMANALO69   | VOLCANO69     | 0.00984 | 0.02975 |  47.9 | —  |
| 56–59  | 19 | 23 | PEARL CITY69  | WAIPAHU69     | 0.01258 | 0.03045 | 100.3 | —  |
| 60     | 20 | 23 | MILILANI69    | WAIPAHU69     | 0.04266 | 0.07484 |  85.8 | —  |
| 61     | 20 | 30 | MILILANI69    | WAHIAWA69     | 0.02780 | 0.06849 |  64.4 | —  |
| 62     | 20 | 32 | MILILANI69    | HAUULA69      | 0.01494 | 0.03992 |  52.0 | —  |
| 63–65  | 22 | 23 | WAIPAHU138    | WAIPAHU69     | 0.00445 | 0.11658 | 149.6 | 1.0 (xfmr) |
| 66     | 22 | 25 | WAIPAHU138    | EWA BEACH138  | 0.00298 | 0.01768 | 277.0 | —  |
| 67     | 22 | 37 | WAIPAHU138    | KAHE138       | 0.00682 | 0.04084 | 292.9 | —  |
| 68–69  | 23 | 34 | WAIPAHU69     | SCHOFIELD69   | 0.03014 | 0.09172 |  74.4 | —  |
| 70–71  | 24 | 26 | KAPOLEI69     | EWA BEACH69   | 0.01645 | 0.05348 |  65.6 | —  |
| 72     | 24 | 33 | KAPOLEI69     | WAIANAE69     | 0.01690 | 0.03319 |  69.5 | —  |
| 73     | 24 | 36 | KAPOLEI69     | COGEN69       | 0.01893 | 0.03755 |  48.7 | —  |
| 74–75  | 25 | 26 | EWA BEACH138  | EWA BEACH69   | 0.00532 | 0.14727 | 131.7 | 1.0 (xfmr) |
| 76–77  | 25 | 35 | EWA BEACH138  | KALAELOA138   | 0.00413 | 0.02442 | 201.4 | —  |
| 78–79  | 25 | 37 | EWA BEACH138  | KAHE138       | 0.00489 | 0.02615 | 224.0 | —  |
| 80     | 27 | 29 | KAHUKU69      | LAIE69        | 0.02388 | 0.04283 |  67.5 | —  |
| 81–82  | 28 | 30 | HALEIWA69     | WAHIAWA69     | 0.01120 | 0.03140 |  70.0 | —  |
| 83     | 28 | 31 | HALEIWA69     | WAIALUA69     | 0.05107 | 0.10468 |  58.3 | —  |
| 84     | 29 | 30 | LAIE69        | WAHIAWA69     | 0.01479 | 0.03666 |  49.3 | —  |
| 85     | 29 | 32 | LAIE69        | HAUULA69      | 0.02315 | 0.06135 |  60.0 | —  |
| 86     | 31 | 33 | WAIALUA69     | WAIANAE69     | 0.02696 | 0.05761 |  38.0 | —  |
| 87–88  | 33 | 34 | WAIANAE69     | SCHOFIELD69   | 0.03528 | 0.08393 |  84.3 | —  |
| 89     | 35 | 37 | KALAELOA138   | KAHE138       | 0.00152 | 0.00837 | 297.6 | —  |

---

# Hawaii case hawaii.json file
> examples/PhasorDynamics/Medium/Hawaii/hawaii.json

## Buses (37 buses)

37 buses numbered 1–37, with identical names and base voltages as in the `.m` file:

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"number"` | `BUS_I` (col 1) | ✓ 1–37, identical |
| `"name"` (e.g. `"ALOHA138"`) | `mpc.bus_name` | ✓ identical |
| `"v_base"` (138000 or 69000) | `BASE_KV` (col 10) | ✓ identical |
| `"init": {Vr, Vi}` | `VM`/`VA` cols (powerflow solution) | ✓ matches powerflow initial conditions |
| `"mon": ["Vr", "Vi"]` | — (GridKit-specific monitoring config) | n/a |
| *(absent)* | `BUS_TYPE` (col 2, PQ/PV/ref) | ✗ not in JSON |
| *(absent)* | `PD`/`QD` (cols 3–4, load) | ✗ not in JSON — loads stored separately as `"class": "Load"` devices |
| *(absent)* | `GS`/`BS` (cols 5–6, shunt conductance/susceptance) | ✗ not in JSON — shunts stored as `"class": "Load"` with X≠0 |
| *(absent)* | `AREA`, `ZONE` (cols 7, 11) | ✗ not in JSON |
| *(absent)* | `VMAX`/`VMIN` (cols 12–13, voltage limits) | ✗ not in JSON |

## Branches (89 branches)

89 `Branch` devices in the `"devices"` array, IDs of the form `BR_<from>_<to>_<parallel_index>` (e.g. `BR_1_2_1`, `BR_1_2_2`, `BR_1_2_3` for the 3 parallel 1→2 circuits). Full 1:1 correspondence with `mpc.branch` — same 89 rows, same R/X/B values.

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"id"` (e.g. `"BR_1_2_1"`) | row index + `F_BUS`/`T_BUS` + parallel count | ✓ derived, not stored directly |
| `"ports": {"bus1", "bus2"}` | `F_BUS` (col 1), `T_BUS` (col 2) | ✓ identical |
| `"params": {"R"}` | `BR_R` (col 3) | ✓ identical |
| `"params": {"X"}` | `BR_X` (col 4) | ✓ identical |
| `"params": {"B"}` | `BR_B` (col 5, total charging susceptance) | ✓ identical |
| `"params": {"G"}` | — (always 0.0 in this case) | ✓ (shunt conductance, not in standard MATPOWER) |
| *(absent)* | `RATE_A/B/C` (cols 6–8, thermal ratings) | ✗ not in JSON |
| *(absent)* | `TAP` (col 9, transformer tap ratio) | ✗ not in JSON — see note below |
| *(absent)* | `SHIFT` (col 10, phase shift) | ✗ not in JSON |
| *(absent)* | `BR_STATUS` (col 11) | ✗ not in JSON (all branches assumed in-service) |
| *(absent)* | Powerflow results (PF, QF, PT, QT cols) | ✗ not in JSON |

**Note on transformers:** The `.m` file marks transformers with a non-zero `TAP` (e.g. `1.00000`). This information is **not preserved** in the JSON — transformers and lines are both stored as `"class": "Branch"` with only R/X/B params. The 13 transformer branches (e.g. 1→2, 5→6, 9→10, 14→15, 22→23, 25→26) lose their tap ratio distinction in the JSON.

## Genrou generators (39 of 45 from .m)

**Total Genrou devices:** 39  
**Unique buses with generators:** 10

### Discrepancy vs. .m file

The `.m` file has **45 generators** in `mpc.gen`; `hawaii.json` has only **39 Genrou** devices. The 6 missing generators were intentionally excluded — `hawaii.json` states: *"A modified Hawaii 37 Bus Case that only has synchronous machines."*

The 6 excluded generators all have `GEN_STATUS = 0` (offline, `PG = 0`) in the `.m` file, and are large-MVA units consistent with inverter-based resources (wind/solar) that have no synchronous machine dynamics model.

### Field mapping: Genrou params vs. .m file

| Field in JSON | Corresponds to `.m` | Match |
|---|---|---|
| `"id"` (e.g. `"2_1"`) | `GEN_BUS` + rank within bus group in `mpc.gen` | ✓ derived; rank reflects row order in `mpc.gen`, gaps where offline gens were removed |
| `"ports": {"bus"}` | `GEN_BUS` (col 1 of `mpc.gen`) | ✓ identical |
| `"params": {"p0"}` | `PG` (col 2) ÷ `mva` (machine MVA base) | ✓ per-unit on machine base |
| `"params": {"q0"}` | `QG` (col 3) ÷ `mva` | ✓ per-unit on machine base |
| `"params": {"mva"}` | `MBASE` (col 7, machine MVA base) | ✓ identical |
| `"params": {"H", "D", "Ra", "Xd", …}` | `.dyr` dynamics file (not in `.m`) | ✓ dynamics params from PowerWorld `.dyr` export |
| `"mon": ["delta", "omega"]` | — (GridKit-specific monitoring config) | n/a |
| *(absent)* | `QMAX`/`QMIN` (cols 4–5, reactive limits) | ✗ not in JSON |
| *(absent)* | `VG` (col 6, voltage setpoint) | ✗ not in JSON |
| *(absent)* | `PMAX`/`PMIN` (cols 9–10, active limits) | ✗ not in JSON |
| *(absent)* | `GEN_STATUS` (col 8) | ✗ not in JSON — offline gens are simply absent |

### Excluded generators

| gen_row | Bus | Bus name    | PMAX (MW) | JSON id (absent) |
|--------:|----:|:------------|----------:|:-----------------|
| 30 | 35 | KALAELOA138 | 180.0 | `35_3` |
| 33 | 35 | KALAELOA138 |  85.0 | `35_6` |
| 40 | 37 | KAHE138     | 128.7 | `37_1` |
| 41 | 37 | KAHE138     | 128.1 | `37_2` |
| 43 | 37 | KAHE138     |  82.1 | `37_4` |
| 45 | 37 | KAHE138     |  77.9 | `37_6` |

This is why bus 35 Genrou IDs skip `_3` and `_6`, and bus 37 has only `37_3` and `37_5` — the rank suffixes reflect original row position within each bus group in `mpc.gen`, not a renumbered compact sequence.

Bus numbers are model indices (1-based), **not** voltage levels. The `138`/`69` in bus names like `ALOHA138` refer to kV voltage level, not the bus number.

| Bus # | Bus name | # gens | Generator IDs | H values | Hops from fault (ALOHA138) |
|-------|----------|--------|---------------|----------|----------------------------|
| 2   | ALOHA69      | 4  | `2_1`, `2_2`, `2_3`, `2_4` | 3.69 (×4) | 1 (transformer to ALOHA138, X=0.124) |
| 23  | WAIPAHU69    | 10 | `23_1`–`23_9`, `23_10` | 6.15 (×8), 3.0 (×2) | 2 (via WAIPAHU138, X≈0.024+0.117) |
| 26  | EWA BEACH69  | 2  | `26_1`, `26_2` | 3.0 (×2) | 2 (via EWA BEACH138, X≈0.016+0.147) |
| 27  | KAHUKU69     | 2  | `27_1`, `27_2` | 3.0 (×2) | far (north shore) |
| 28  | HALEIWA69    | 2  | `28_1`, `28_2` | 3.0 (×2) | far (north shore) |
| 33  | WAIANAE69    | 1  | `33_1` | 3.0 | 3+ (west side) |
| 34  | SCHOFIELD69  | 6  | `34_1`–`34_6` | 4.35 (×6) | 3+ (central Oahu) |
| 35  | KALAELOA138  | 6  | `35_1`, `35_2`, `35_4`, `35_5`, `35_7`, `35_8` | 5.22 (×6) | 1 (direct branch, X=0.042) |
| 36  | COGEN69      | 4  | `36_1`–`36_4` | 4.42 (×4) | 3+ |
| 37  | KAHE138      | 2  | `37_3`, `37_5` | 2.48 (×2) | 2 (via KALAELOA138, lowest path X≈0.050) |

## Device model triplet: Genrou + Tgov1 + Ieeet1

Every generator in the Hawaii case is modeled as a triplet of three coupled devices:

- **Genrou** — the synchronous machine itself (rotor electromechanics). Tracks rotor angle `delta`, rotor speed `omega`, and internal flux linkages. Inertia `H` is the epistemic uncertainty parameter used in the current UQ campaign; other Genrou parameters (e.g. `D`, `Xd`, `Xdp`) are candidates for future epistemic UQ runs. Aleatoric uncertainty sources (e.g. load variability, renewable output) will be added in later campaigns.
- **Tgov1** — a simplified steam/gas turbine governor. Receives a speed error signal and outputs mechanical power `pmech` back to Genrou. Without it, mechanical power would be fixed, and the machine would not respond to frequency deviations.
- **Ieeet1** — an IEEE Type 1 excitation system / automatic voltage regulator (AVR). Senses terminal voltage and regulates field voltage `efd` into Genrou. Without it, field voltage would be fixed, and the machine would not regulate voltage.

The three are wired together via named ports: Genrou exposes `pmech` and `efd` input ports; Tgov1 and Ieeet1 write to those ports respectively. All 39 Genrou devices have a complete paired Tgov1 and Ieeet1 — no generators are modeled as infinite-bus or constant-parameter machines.

### Port wiring and grid connection

Each device has a `ports` block with two distinct kinds of entries:

| Port key | Type | Value | Meaning |
|----------|------|-------|---------|
| `"bus"` | electrical | bus number (int) | connects the machine to that network bus — this is the only grid connection |
| `"speed"` | signal | global index (int) | rotor speed signal, shared between Genrou and its Tgov1 and Ieeet1 |
| `"pmech"` | signal | global index (int) | mechanical power signal, shared between Genrou and its Tgov1 |
| `"efd"` | signal | global index (int) | field voltage signal, shared between Genrou and its Ieeet1 |

Signal port indices are **globally unique across all devices** in the file. For each generator triplet, the three devices share the same `speed`, `pmech`, `efd` index values — that is how they are wired together. The `bus` port value is just the bus number; the bus itself lives in the `buses` array and is never modified when removing a generator.

**Bus 2 has 4 generators.** Their signal port allocations show the pattern clearly:

| id | class | bus | speed | pmech | efd |
|----|-------|-----|-------|-------|-----|
| `2_1` | Genrou | 2 | 0 | 1 | 78 |
| `2_1_tgov1` | Tgov1 | 2 | 0 | 1 | — |
| `2_1_ieeet1` | Ieeet1 | 2 | 0 | — | 78 |
| `2_2` | Genrou | 2 | 2 | 3 | 79 |
| `2_2_tgov1` | Tgov1 | 2 | 2 | 3 | — |
| `2_2_ieeet1` | Ieeet1 | 2 | 2 | — | 79 |
| `2_3` | Genrou | 2 | 4 | 5 | 80 |
| `2_3_tgov1` | Tgov1 | 2 | 4 | 5 | — |
| `2_3_ieeet1` | Ieeet1 | 2 | 4 | — | 80 |
| `2_4` | Genrou | 2 | 6 | 7 | 81 |
| `2_4_tgov1` | Tgov1 | 2 | 6 | 7 | — |
| `2_4_ieeet1` | Ieeet1 | 2 | 6 | — | 81 |

**Removing a generator** means removing all 3 devices (Genrou + Tgov1 + Ieeet1) from the `devices` array. The bus node itself is untouched. No separate "connection" device is involved. The signal port indices freed up by removal do not need to be re-numbered — GridKit resolves ports by matching index values within the device list, so gaps are fine as long as no other device references the removed indices.

### Example JSON components (bus 2, id `2_1`)
* machine Genrou
* governor Tgov1
* excitor Ieeet1


Below are trimmed excerpts from `hawaii.json` showing the key:value structure, especially the `params` blocks.

```json
{
	"class": "Genrou",
	"ports": {
		"bus": 2,
		"speed": 0,
		"pmech": 1,
		"efd": 78
	},
	"id": "2_1",
	"params": {
		"p0": 0.025000000400000003,
		"q0": 0.0080000004,
		"H": 3.69000006,
		"D": 0.0,
		"Ra": 0.0,
		"Tdop": 7.61999989,
		"Tdopp": 0.06,
		"Tqop": 1.30999994,
		"Tqopp": 0.04,
		"Xd": 1.54999995,
		"Xdp": 0.15000001,
		"Xdpp": 0.15000001,
		"Xq": 1.51999998,
		"Xqp": 0.47999999,
		"Xqpp": 0.15000001,
		"Xl": 0.12,
		"S10": 0.17,
		"S12": 0.56999999,
		"mva": 2.79999995
	},
	"mon": ["delta", "omega"]
}
```

```json
{
	"class": "Tgov1",
	"ports": {
		"bus": 2,
		"speed": 0,
		"pmech": 1
	},
	"id": "2_1_tgov1",
	"params": {
		"R": 0.05,
		"T1": 0.5,
		"T2": 2.5,
		"T3": 7.5,
		"Pvmin": 0.0,
		"Pvmax": 1.0,
		"Dt": 0.0
	}
}
```

```json
{
	"class": "Ieeet1",
	"ports": {
		"bus": 2,
		"speed": 0,
		"efd": 78
	},
	"id": "2_1_ieeet1",
	"params": {
		"Tr": 0.01,
		"Ka": 24.79,
		"Ta": 0.02,
		"Ke": 1.0,
		"Te": 0.72,
		"Kf": 0.09,
		"Tf": 1.2,
		"Vrmin": -100.0,
		"Vrmax": 100.0,
		"E1": 2.8,
		"E2": 3.73,
		"Se1": 0.04,
		"Se2": 0.33,
		"Ispdlim": 0.0
	}
}
```

## UQ parameter selection 

One generator per bus, 4 buses chosen for diverse H values:

| id | bus # | bus name | H (nominal) |
|----|-------|----------|-------------|
| `2_1`  | 2  | ALOHA69     | 3.69 |
| `23_1` | 23 | WAIPAHU69   | 6.15 |
| `34_1` | 34 | SCHOFIELD69 | 4.35 |
| `35_1` | 35 | KALAELOA138 | 5.22 |

### v1 (20 samples)

- dist: uniform ±10% of nominal
- method: LHS, seed=42
- serialize: stacked (`results.parquet`)

### v2 (1K samples)

- dist: uniform ±10% of nominal
- method: LHS, seed=42
- serialize: per-run (`runs/run_NNN.parquet`)

| id | H nominal | lo | hi |
|----|-----------|----|----|
| `2_1`  | 3.69 | 3.321 | 4.059 |
| `23_1` | 6.15 | 5.535 | 6.765 |
| `34_1` | 4.35 | 3.915 | 4.785 |
| `35_1` | 5.22 | 4.698 | 5.742 |

### v3 (1K samples) — epistemic UQ

- dist: **Gaussian**, std = 12% of nominal (epistemic uncertainty in H)
- method: LHS + normal ppf (scipy), seed=42
- serialize: per-run (`runs/run_NNN.parquet`)
- data path: `/kfs2/projects/scidac/scidac-data/gridkit-runs/hawaii-v3/`
- metadata: `hawaii-v3/meta.yml`, samples: `hawaii-v3/samples.csv`

| id | H nominal | mean | std (12%) |
|----|-----------|------|-----------|
| `2_1`  | 3.69 | 3.69 | 0.4428 |
| `23_1` | 6.15 | 6.15 | 0.7380 |
| `34_1` | 4.35 | 4.35 | 0.5220 |
| `35_1` | 5.22 | 5.22 | 0.6264 |

### v4 (4K samples) — epistemic UQ, 2 generators

- dist: **Gaussian**, std = 12% of nominal (epistemic uncertainty in H)
- method: LHS + normal ppf (scipy), seed=42
- serialize: per-run (`runs/run_NNN.parquet`)
- generators: `2_1` and `23_1` only (`34_1`, `35_1` excluded)
- data path: `/kfs2/projects/scidac/scidac-data/gridkit-runs/hawaii-v4/`
- metadata: `hawaii-v4/meta.yml`, samples: `hawaii-v4/samples.csv`

| id | H nominal | mean | std (12%) |
|----|-----------|------|-----------|
| `2_1`  | 3.69 | 3.69 | 0.4428 |
| `23_1` | 6.15 | 6.15 | 0.7380 |

### v5 (16K samples) — epistemic UQ, 4 generators

- dist: **Gaussian**, std = 12% of nominal (epistemic uncertainty in H)
- method: LHS + normal ppf (scipy), seed=42
- serialize: per-run (`runs/run_NNN.parquet`)
- generators: all 4 (`2_1`, `23_1`, `34_1`, `35_1`)
- solver overrides: none (base case hawaii.solver.json)
- monitors: `Vm`, `Va` (buses, polar form of `Vr`/`Vi`); `delta`, `omega` (genrou)
- data path: `/kfs2/projects/scidac/scidac-data/gridkit-runs/hawaii-v5/`
- metadata: `hawaii-v5/meta.yml`, samples: `hawaii-v5/samples.csv`

| id | bus | bus name | H nominal | mean | std (12%) |
|----|-----|----------|-----------|------|-----------|
| `2_1`  | 2  | ALOHA69     | 3.69 | 3.69 | 0.4428 |
| `23_1` | 23 | WAIPAHU69   | 6.15 | 6.15 | 0.7380 |
| `34_1` | 34 | SCHOFIELD69 | 4.35 | 4.35 | 0.5220 |
| `35_1` | 35 | KALAELOA138 | 5.22 | 5.22 | 0.6264 |


