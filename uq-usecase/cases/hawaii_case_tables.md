# Hawaii case JSON tables (hawaii.json)

> `examples/PhasorDynamics/Medium/Hawaii/hawaii.json`

Parsed directly from the JSON file. All values are exact as stored in the file.
Powers (`p0`, `q0`) are in per-unit on each machine's own MVA base.
Reactances (`Xd`, `Xdp`, etc.) are in per-unit on machine base.
`H` is in seconds. `mva` is machine MVA base.
`H*mva` (MWs) is the total stored kinetic energy — the physically meaningful inertia contribution to the grid.

## Genrou devices (39 total)

Sorted by `H*mva` ascending. Bus names are as in `hawaii.json` buses array.

**Why H×mva?** `H` (the inertia constant, in seconds) is defined as stored kinetic energy normalized by the machine's own MVA rating: $H = \frac{\frac{1}{2} J \omega_0^2}{S_\text{rated}}$. It is already per-unit on each machine's individual base, not on the system base (100 MVA). Each machine has its own `mva` (`MBASE` in MATPOWER), which varies widely across the fleet — from 2.8 MVA for the small ALOHA69 units to 124.3 MVA for the largest KALAELOA138 unit. To recover the actual stored kinetic energy in MWs you multiply back: $H \times \text{mva}$. Two machines with the same H can contribute very different inertia to the grid: `2_1` (H=3.69, mva=2.8) stores 10.3 MWs, while `23_1` (H=6.15, mva=85.6) stores 526 MWs. **H×mva is therefore the physically relevant quantity for ranking inertia contributions and selecting UQ generators.**

### Genrou parameter glossary

All quantities are in per-unit on the machine's own MVA base unless noted otherwise.

| param | unit | description |
|-------|------|-------------|
| `mva` | MVA | machine MVA rating (each machine's own base; varies per unit) |
| `p0` | pu | initial active power output (per-unit on `mva`) |
| `q0` | pu | initial reactive power output (per-unit on `mva`) |
| `H` | s | inertia constant: stored kinetic energy at rated speed divided by `mva`; epistemic UQ parameter |
| `D` | pu | damping coefficient; zero for all machines in this case |
| `Ra` | pu | armature (stator) resistance; zero for all machines in this case |
| `Xd` | pu | d-axis synchronous reactance |
| `Xdp` | pu | d-axis transient reactance (primed) |
| `Xdpp` | pu | d-axis subtransient reactance (double-primed) |
| `Xq` | pu | q-axis synchronous reactance |
| `Xqp` | pu | q-axis transient reactance |
| `Xqpp` | pu | q-axis subtransient reactance |
| `Xl` | pu | leakage (armature) reactance |
| `Tdop` | s | d-axis transient open-circuit time constant |
| `Tdopp` | s | d-axis subtransient open-circuit time constant |
| `Tqop` | s | q-axis transient open-circuit time constant |
| `Tqopp` | s | q-axis subtransient open-circuit time constant |
| `S10` | pu | saturation factor at 1.0 pu flux linkage |
| `S12` | pu | saturation factor at 1.2 pu flux linkage; together with `S10` defines the saturation curve |

The suffix convention `p` = transient (one rotor winding), `pp` = subtransient (two rotor windings). Open-circuit time constants (`Tdo`, `Tqo`) are for the rotor circuit with stator open; they are longer than short-circuit time constants by the ratio $X / X'$.

| id | bus | bus_name | mva | H | H&times;mva (MWs) | p0 (pu) | q0 (pu) | D | Ra | Xd | Xdp | Xdpp | Xq | Xqp | Xqpp | Xl | Tdop | Tdopp | Tqop | Tqopp | S10 | S12 |
|----|-----|----------|----:|------:|-------:|--------:|--------:|----:|----:|------:|------:|------:|------:|------:|------:|-----:|-----:|------:|-----:|------:|-----:|-----:|
| `2_1` | 2 | ALOHA69 | 2.800 | 3.690 | **10.33** | 0.02500 | 0.00800 | 0.0 | 0.0 | 1.550 | 0.150 | 0.150 | 1.520 | 0.480 | 0.150 | 0.12 | 7.620 | 0.060 | 1.310 | 0.040 | 0.170 | 0.570 |
| `2_2` | 2 | ALOHA69 | 2.800 | 3.690 | **10.33** | 0.02500 | 0.00800 | 0.0 | 0.0 | 1.550 | 0.150 | 0.150 | 1.520 | 0.480 | 0.150 | 0.12 | 7.620 | 0.060 | 1.310 | 0.040 | 0.170 | 0.570 |
| `2_3` | 2 | ALOHA69 | 2.800 | 3.690 | **10.33** | 0.02500 | 0.00800 | 0.0 | 0.0 | 1.550 | 0.150 | 0.150 | 1.520 | 0.480 | 0.150 | 0.12 | 7.620 | 0.060 | 1.310 | 0.040 | 0.170 | 0.570 |
| `2_4` | 2 | ALOHA69 | 2.800 | 3.690 | **10.33** | 0.02500 | 0.00800 | 0.0 | 0.0 | 1.550 | 0.150 | 0.150 | 1.520 | 0.480 | 0.150 | 0.12 | 7.620 | 0.060 | 1.310 | 0.040 | 0.170 | 0.570 |
| `36_1` | 36 | COGEN69 | 3.500 | 4.420 | **15.47** | 0.02000 | 0.01000 | 0.0 | 0.0 | 1.870 | 0.220 | 0.170 | 1.780 | 0.550 | 0.170 | 0.13 | 5.220 | 0.040 | 0.550 | 0.070 | 0.060 | 0.360 |
| `36_2` | 36 | COGEN69 | 3.500 | 4.420 | **15.47** | 0.02101 | 0.01000 | 0.0 | 0.0 | 1.870 | 0.220 | 0.170 | 1.780 | 0.550 | 0.170 | 0.13 | 5.220 | 0.040 | 0.550 | 0.070 | 0.060 | 0.360 |
| `36_3` | 36 | COGEN69 | 3.500 | 4.420 | **15.47** | 0.02000 | 0.01000 | 0.0 | 0.0 | 1.870 | 0.220 | 0.170 | 1.780 | 0.550 | 0.170 | 0.13 | 5.220 | 0.040 | 0.550 | 0.070 | 0.060 | 0.360 |
| `36_4` | 36 | COGEN69 | 3.500 | 4.420 | **15.47** | 0.02480 | 0.01000 | 0.0 | 0.0 | 1.870 | 0.220 | 0.170 | 1.780 | 0.550 | 0.170 | 0.13 | 5.220 | 0.040 | 0.550 | 0.070 | 0.060 | 0.360 |
| `26_2` | 26 | EWA BEACH69 | 11.200 | 3.000 | **33.60** | 0.10200 | 0.03100 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `34_1` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.07220 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `34_2` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.06040 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `34_3` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.06040 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `34_4` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.06040 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `34_5` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.06009 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `34_6` | 34 | SCHOFIELD69 | 9.200 | 4.350 | **40.02** | 0.07220 | -0.00982 | 0.0 | 0.0 | 1.740 | 0.160 | 0.160 | 1.670 | 0.290 | 0.160 | 0.12 | 4.990 | 0.020 | 0.890 | 0.040 | 0.170 | 0.580 |
| `23_9` | 23 | WAIPAHU69 | 16.200 | 3.000 | **48.60** | 0.12532 | 0.00044 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `26_1` | 26 | EWA BEACH69 | 22.000 | 3.000 | **66.00** | 0.20000 | 0.06000 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `27_2` | 27 | KAHUKU69 | 30.400 | 3.000 | **91.20** | 0.27600 | -0.05000 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `33_1` | 33 | WAIANAE69 | 30.400 | 3.000 | **91.20** | 0.27600 | 0.08300 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `27_1` | 27 | KAHUKU69 | 33.000 | 3.000 | **99.00** | 0.30000 | -0.05400 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `23_10` | 23 | WAIPAHU69 | 50.500 | 3.000 | **151.50** | 0.39130 | 0.00044 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `28_1` | 28 | HALEIWA69 | 53.900 | 3.000 | **161.70** | 0.46300 | -0.05487 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `37_5` | 37 | KAHE138 | 85.900 | 2.480 | **213.03** | 0.68480 | 0.05626 | 0.0 | 0.0 | 1.570 | 0.180 | 0.130 | 1.540 | 0.270 | 0.130 | 0.10 | 5.900 | 0.040 | 0.580 | 0.070 | 0.130 | 0.510 |
| `28_2` | 28 | HALEIWA69 | 75.900 | 3.000 | **227.70** | 0.69000 | -0.05487 | 0.0 | 0.0 | 2.100 | 0.200 | 0.180 | 0.500 | 0.500 | 0.180 | 0.15 | 7.000 | 0.040 | 0.750 | 0.050 | 0.000 | 0.000 |
| `37_3` | 37 | KAHE138 | 95.900 | 2.480 | **237.83** | 0.58324 | 0.05626 | 0.0 | 0.0 | 1.570 | 0.180 | 0.130 | 1.540 | 0.270 | 0.130 | 0.10 | 5.900 | 0.040 | 0.580 | 0.070 | 0.130 | 0.510 |
| `35_4` | 35 | KALAELOA138 | 22.000 | 5.220 | **114.84** | 0.20000 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |
| `35_2` | 35 | KALAELOA138 | 30.800 | 5.220 | **160.78** | 0.28000 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |
| `35_7` | 35 | KALAELOA138 | 55.000 | 5.220 | **287.10** | 0.50000 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |
| `23_8` | 23 | WAIPAHU69 | 51.900 | 6.150 | **319.19** | 0.43188 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `23_7` | 23 | WAIPAHU69 | 52.500 | 6.150 | **322.88** | 0.43614 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `35_1` | 35 | KALAELOA138 | 63.800 | 5.220 | **332.84** | 0.58000 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |
| `23_5` | 23 | WAIPAHU69 | 56.300 | 6.150 | **346.25** | 0.44533 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `23_6` | 23 | WAIPAHU69 | 56.300 | 6.150 | **346.25** | 0.44533 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `23_4` | 23 | WAIPAHU69 | 57.000 | 6.150 | **350.55** | 0.47110 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `23_3` | 23 | WAIPAHU69 | 57.100 | 6.150 | **351.17** | 0.47637 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `35_5` | 35 | KALAELOA138 | 93.500 | 5.220 | **488.07** | 0.67000 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |
| `23_1` | 23 | WAIPAHU69 | 85.600 | 6.150 | **526.44** | 0.69275 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `23_2` | 23 | WAIPAHU69 | 85.600 | 6.150 | **526.44** | 0.69275 | 0.00044 | 0.0 | 0.0 | 1.600 | 0.320 | 0.320 | 1.560 | 0.590 | 0.320 | 0.25 | 9.230 | 0.030 | 0.540 | 0.060 | 0.130 | 0.490 |
| `35_8` | 35 | KALAELOA138 | 124.300 | 5.220 | **648.85** | 0.56200 | 0.00053 | 0.0 | 0.0 | 1.790 | 0.400 | 0.280 | 1.710 | 0.670 | 0.280 | 0.22 | 6.240 | 0.030 | 1.060 | 0.080 | 0.150 | 0.530 |

## H and inertia summary by bus

Total stored kinetic energy per bus = sum of H×mva over all Genrou devices at that bus.

| bus | bus_name | n_gens | H (all same?) | mva range (MVA) | H×mva per unit (MWs) | bus total H×mva (MWs) |
|-----|----------|-------:|---------------|-----------------|---------------------:|----------------------:|
| 2 | ALOHA69 | 4 | 3.690 (all same) | 2.800 | 10.33 each | **41.32** |
| 23 | WAIPAHU69 | 10 | 3.000 (×2), 6.150 (×8) | 16.2 – 85.6 | 48.6 – 526.44 | **3,088.66** |
| 26 | EWA BEACH69 | 2 | 3.000 (all same) | 11.2 – 22.0 | 33.6 – 66.0 | **99.60** |
| 27 | KAHUKU69 | 2 | 3.000 (all same) | 30.4 – 33.0 | 91.2 – 99.0 | **190.20** |
| 28 | HALEIWA69 | 2 | 3.000 (all same) | 53.9 – 75.9 | 161.7 – 227.7 | **389.40** |
| 33 | WAIANAE69 | 1 | 3.000 | 30.400 | 91.20 | **91.20** |
| 34 | SCHOFIELD69 | 6 | 4.350 (all same) | 9.200 (all same) | 40.02 each | **240.12** |
| 35 | KALAELOA138 | 6 | 5.220 (all same) | 22.0 – 124.3 | 114.84 – 648.85 | **2,032.48** |
| 36 | COGEN69 | 4 | 4.420 (all same) | 3.500 (all same) | 15.47 each | **61.88** |
| 37 | KAHE138 | 2 | 2.480 (all same) | 85.9 – 95.9 | 213.03 – 237.83 | **450.86** |

**Total system H×mva: 6,685.72 MWs** (sum over all 39 Genrou devices, computed from table above)

## H values across the fleet

There are only 4 distinct H values in this case:

| H (s) | devices | buses |
|------:|--------:|-------|
| 2.480 | 2 | 37 (KAHE138) |
| 3.000 | 13 | 23 (×2), 26 (×2), 27 (×2), 28 (×2), 33 (×1) |
| 3.690 | 4 | 2 (ALOHA69) |
| 4.350 | 6 | 34 (SCHOFIELD69) |
| 4.420 | 4 | 36 (COGEN69) |
| 5.220 | 6 | 35 (KALAELOA138) |
| 6.150 | 8 | 23 (WAIPAHU69) |

> **Note on "smallest H" vs "smallest inertia":** These are not the same selection.
> - 4 smallest H: `37_3`, `37_5` (H=2.48), then two of the 13 H=3.0 gens (ambiguous — many candidates).
> - 4 smallest H×mva: `2_1`, `2_2`, `2_3`, `2_4` (all at ALOHA69, H×mva=10.33 MWs each) — unambiguous but all 4 are identical machines at the same bus, which limits UQ diversity.
> - For UQ gen selection, "one representative per bus, smallest H×mva" is likely more useful than raw ranking.

## Distinct parameter groups

Within each bus, all Genrou at the same bus share identical dynamics params (H, D, Ra, Xd, Xdp, etc.) — only `p0`, `q0`, and `mva` differ slightly between units. There are 6 distinct parameter groups (one per bus):

| group | buses | H | Xd | Xdp | Xq | Xqp | Tdop | Tqop |
|-------|-------|---:|----|-----|-----|-----|------|------|
| ALOHA69 type | 2 | 3.690 | 1.550 | 0.150 | 1.520 | 0.480 | 7.620 | 1.310 |
| WAIPAHU69 type (H=6.15) | 23 (×8) | 6.150 | 1.600 | 0.320 | 1.560 | 0.590 | 9.230 | 0.540 |
| H=3.0 type | 23 (×2), 26, 27, 28, 33 | 3.000 | 2.100 | 0.200 | 0.500 | 0.500 | 7.000 | 0.750 |
| SCHOFIELD69 type | 34 | 4.350 | 1.740 | 0.160 | 1.670 | 0.290 | 4.990 | 0.890 |
| KALAELOA138 type | 35 | 5.220 | 1.790 | 0.400 | 1.710 | 0.670 | 6.240 | 1.060 |
| COGEN69 type | 36 | 4.420 | 1.870 | 0.220 | 1.780 | 0.550 | 5.220 | 0.550 |
| KAHE138 type | 37 | 2.480 | 1.570 | 0.180 | 1.540 | 0.270 | 5.900 | 0.580 |

## BusFault device

The case has one `BusFault` device, at bus 1 (ALOHA138), active at $t=0$ (`state0: true`).
This is distinct from the Illinois case where faults are at every bus but disabled by default.

## Solver (hawaii.solver.json)

Base solver configuration (before any UQ overrides):

| parameter | value |
|-----------|-------|
| tmax | [read from hawaii.solver.json] |
| fault bus | 1 (ALOHA138) |
| fault on | [read from hawaii.solver.json] |
| fault off | [read from hawaii.solver.json] |

> Values marked `[read from ...]` are not stored here to avoid staleness; read the file directly or check `meta.yml` from any run.
