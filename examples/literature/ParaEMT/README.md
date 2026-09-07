# ParaEMT / GridKit 9-bus comparison

The full governor-step and ideal generator-trip cases run in GridKit.
The required SEXS-PTI, GASTPTI, and IEEEST models were ported from
PhasorDynamics and committed separately; all 20 EMT tests pass.
Cases, executable driver, simulation outputs, and reproducible plots are here.

**[Open the plot gallery](plots/index.html)** ·
[Trip PDF](plots/trip/comparison.pdf) ·
[Governor-step PDF](plots/governor_step/comparison.pdf) ·
[Runtime comparison](RUNTIMES.md)

Each figure contains one response: **ParaEMT on top, GridKit in the middle,
GridKit − ParaEMT on the bottom**, with identical time and vertical limits
in all three panels. Three-phase signals are plotted together. Bus waveform
views show 0.990–1.025 s, approximately two cycles around the event, plus a
separate late-time view. Slow responses retain full-run and event views.
Speed is displayed as deviation from synchronous speed, in pu.

There are **60 figures per experiment**, covering all 48 normalized channels:
27 phase voltages, nine voltage magnitudes, and speed, mechanical power,
field voltage, and stabilizer output for each of the three generators.
`plots/coverage.json` records every figure's channels, shared axes and limits.
Earlier plot filenames now point to figures in this layout.

- [Bus 1 three-phase trip waveform](plots/trip/bus1_abc_event.png)
- [Bus 4 three-phase trip waveform](plots/trip/bus4_abc_event.png)
- [Generator 2 trip speed response](plots/trip/gen2_speed_pu_full.png)
- [Generator 1 trip speed response](plots/trip/gen1_speed_pu_full.png)
- [Generator 1 governor-step speed response](plots/governor_step/gen1_speed_pu_full.png)

## Experiments and interpretation

Both events occur at 1 s; each run covers 0–3 s. Every generator retains
its full GridKit machine, governor, exciter and stabilizer. Each simulator
saves 60,001 common samples spaced 50 µs apart. ParaEMT also saves every
integration step between 0.995 and 1.015 s in `event_waveforms.csv.gz`.
The bus event plots use these native samples inside that interval and
regular samples outside it. Error panels always use the common 50 µs
grid, without interpolation, phase alignment or fitted gains.

| Experiment | GridKit | ParaEMT |
| --- | --- | --- |
| Governor step | G1 governor reference −0.02 pu on its machine base | Same supplied GAST reference step |
| Generator trip | Ideal opening between G1 and bus 1, with a separate isolated machine terminal | Upstream `GenTrip` removes G1 Norton injection and conductance; its electrical history stops updating while its state kernel continues |

**The post-trip G1 models differ.** GridKit's disconnected generator
accelerates and its governor responds; ParaEMT's archived G1 internal
states are not an equivalent isolated-machine trajectory. All G1 traces
remain visible. Mechanical power is the turbine output, not electrical
power injected into the network after disconnection.

The GridKit opening explicitly projects stator and transformer currents
to zero while preserving rotor winding fluxes, rotor motion, controller
states, and all other network differential states. It represents finite
pre/post-event states, without an arc, snubber, or finite voltage-impulse
amplitude. The common-grid GridKit sample at 1 s is the **left limit**;
ParaEMT's sample includes its trip step. Both GridKit monitor limits and
the exact state projection are archived.

ParaEMT's bus-1 trip voltage magnitude at 1 s grows from **3.595605 pu**
to **8.062099 pu** to **17.045898 pu** at 50, 25 and 12.5 µs. The spike is
retained in all relevant plots; it is not a converged transient-voltage
answer key. These trip comparisons demonstrate responses and numerical
differences, not complete validation of an identical breaker model.
[GRIDKIT_EQUIVALENCE.md](GRIDKIT_EQUIVALENCE.md) gives the mapping and event equations.

## Measured agreement

Maximum absolute governor-step discrepancies over all saved samples and
channels in each group, using GridKit tolerance `1e-9` and `mu=50000`:

| Quantity, pu | vs ParaEMT 50 µs | vs 25 µs | vs 12.5 µs | GridKit 1e-8 vs 1e-9 |
| --- | ---: | ---: | ---: | ---: |
| Phase voltage | 0.0256288 | 0.0129372 | 0.00653278 | 5.34048e-07 |
| Voltage magnitude | 0.000351683 | 0.000176497 | 0.000303192 | 1.28087e-07 |
| Rotor speed | 3.01155e-05 | 1.51938e-05 | 7.66617e-06 | 1.01776e-09 |
| Mechanical power | 0.000322696 | 0.000162724 | 8.20427e-05 | 7.0769e-09 |
| Field voltage | 0.000392251 | 0.000196044 | 9.7868e-05 | 2.1424e-08 |
| Stabilizer output | 0.00132954 | 0.000677736 | 0.000349688 | 1.31205e-06 |

Most differences approximately halve under ParaEMT time-step refinement.
The finest voltage-magnitude maximum occurs at bus 3 at 50 µs during
startup. The original power flow requires a voltage-phasor adjustment
of at most `4.85e-6 pu` to give a consistent continuous inductive initial
state. Initialization and numerical-method differences remain.

For the trip, excluding G1 internal signals and considering samples after
1.02 s, the finest comparison has maximum differences of `1.15e-5 pu`
in G2/G3 speed, `1.10e-4 pu` in mechanical power and field voltage,
and `3.15e-4 pu` in stabilizer output. Voltage-magnitude discrepancy is
larger (`0.03076 pu`, bus 1 at 1.02385 s). The full intervals, all channels,
all three ParaEMT steps, and both GridKit tolerance refinements are in
[governor metrics](results/gridkit_comparison.json) and
[trip metrics](results/gridkit_trip_comparison.json). There is no
literature-defined pass/fail tolerance for these locally generated data.

## Reproduce

Use Python 3.12 and a built GridKit checkout with Enzyme and SUNDIALS KLU.
The standalone CMake project imports the existing GridKit build's libraries.
From this directory:

```bash
python3 -m venv /tmp/paraemt-reference-venv
/tmp/paraemt-reference-venv/bin/pip install -r upstream/requirements.txt
```

The six saved ParaEMT references use `--event governor-step` or `--event trip`,
with `--dt-us 50`, `25`, or `12.5`. For example:

```bash
/tmp/paraemt-reference-venv/bin/python run_reference.py --event trip --dt-us 12.5 --save-states --output results/trip/dt12_5us > results/trip/dt12_5us/run.log 2>&1
MPLCONFIGDIR=/tmp/paraemt-matplotlib /tmp/paraemt-reference-venv/bin/python run_gridkit.py --gridkit-build ../../../build/emt-rational
python3 verify.py --regenerated
```

`run_gridkit.py` generates both cases, builds the driver, runs each event
at three IDA tolerances, normalizes and compares the data, and regenerates
the figures. To regenerate only figures from saved CSVs:

```bash
MPLCONFIGDIR=/tmp/paraemt-matplotlib /tmp/paraemt-reference-venv/bin/python plot_comparison.py
```

Each GridKit result folder includes compressed normalized and raw monitor
CSVs, a column schema, run log, solver statistics, executable hash,
initial/final assembled-Jacobian checks, both event state limits, and
monitor restart adjustments. The `1e-9` runs also retain full state and
derivative CSVs with explicit column maps. Finest sampled original KCL
mismatches are `0.000162 A` for the governor event and `0.000087 A` for
the trip, independently recomputed by `verify.py`.

`benchmark.py` measures three runs per setting with result capture disabled,
pinned to one logical CPU. ParaEMT kernels are warmed on a disposable state
before timing, with a check that no JIT specialization occurs in the timed
loop. The runner checks each final state against the captured
trajectory. It reports loop wall/CPU times, initialization, process wall time,
and accepted/fixed step counts in [RUNTIMES.md](RUNTIMES.md).
[Cold-process costs](RUNTIMES_COLD.md) are retained separately.

The driver handles initialization, writable setpoints, terminal constraints,
and the ideal opening. These operations are not yet available together in
the ordinary `EMTDynamicSimulation` application.

## Source and provenance

References were generated locally from the authors'
[revision d79d735a](https://github.com/NatLabRockies/ParaEMT_public/tree/d79d735a4a587d56c5b88187d1a499195b6b2b84).
They are not published waveform answer keys or a reproduction of a specific
paper figure. The case contains nine buses, six pi lines, three transformer
series branches, three loads, and three each of GENROU, SEXS, GAST and IEEEST;
no IBRs. It uses serial sparse LU, constant-impedance load option 1, and
bus measurements/PLL enabled from zero.

Original inputs are `ieee9.raw`, `pfd_9_1_1.json`, and `9bus.xlsx`.
Root-level `gen.csv`, `exc.csv`, `gov.csv`, and `pss.csv` are worksheet
exports. Original source is preserved byte for byte under `upstream/`,
with its BSD-3-Clause [license](LICENSE.md). The wrapper preserves the
upstream stepping order, trip history handling, and model equations.
An `ElementTree.getiterator` compatibility alias lets pinned `xlrd==1.2.0`
read XLSX on current Python. The serial path does not require METIS.
[environment.txt](environment.txt) records the environment.

[sources.json](sources.json) records SHA-256 hashes and upstream Git blob
identities, including hashes of local raw outputs. Normalized waveform CSVs,
event samples, solver-work data, run summaries, and plots are tracked. Logs,
raw monitor CSVs, and full-state archives are ignored and retained locally;
a fresh checkout must regenerate them before running the full verification.
`python3 verify.py` checks the frozen artifact snapshot,
complete data, native-step samples, event invariants, independent original
KCL, winding inverses, plot channel coverage, PDF page counts, and gallery
links. Use `--regenerated` after rerunning to skip frozen local hashes
while retaining source and numerical checks.

Reference: Xiong et al., *ParaEMT: An Open Source, Parallelizable, and
HPC-Compatible EMT Simulator for Large-Scale IBR-Rich Power Grids*, IEEE
Transactions on Power Delivery 39(2), 911–921 (2024),
[DOI](https://doi.org/10.1109/TPWRD.2023.3342715).
