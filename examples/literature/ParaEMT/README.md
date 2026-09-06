# ParaEMT 9-bus reference simulations

**GridKit equivalence is blocked by missing EMT SEXS, GAST and IEEEST
controllers. No equivalent GridKit case was created or simulated.** The
corresponding PhasorDynamics classes use different component and signal
interfaces. [GRIDKIT_EQUIVALENCE.md](GRIDKIT_EQUIVALENCE.md) records the
required work, source locations, and machine/network parameter mapping.

Four ParaEMT simulations were completed on 2026-09-06 using the authors'
[revision d79d735a](https://github.com/NatLabRockies/ParaEMT_public/tree/d79d735a4a587d56c5b88187d1a499195b6b2b84).
These are locally generated reference candidates, not published waveform
data or a reproduction of a particular figure in the paper. Model equations
and imported upstream source files are unchanged.

| Experiment | Settings | Results |
| --- | --- | --- |
| Governor step | Generator at bus 1: GAST reference −0.02 pu at 1 s; generators remain connected | [Response plot](plots/governor_step_response.png), [50 µs CSV](results/governor_step/dt50us/reference.csv.gz), [25 µs CSV](results/governor_step/dt25us/reference.csv.gz) |
| Generator trip | Upstream `GenTrip`: generator at bus 1 disconnected at 1 s | [Response plot](plots/trip_response.png), [trip waveform plot](plots/trip_waveforms.png), [50 µs CSV](results/trip/dt50us/reference.csv.gz), [25 µs CSV](results/trip/dt25us/reference.csv.gz) |

Each run covers 0–3 s with serial sparse LU, load option 1 (constant RLC),
and ParaEMT's bus measurements/PLL enabled from time zero. The case contains
9 buses, 6 pi lines, 3 transformer series branches, 3 loads, and 3 each of
GENROU, SEXS, GAST and IEEEST; no IBRs. Initialization uses the supplied
solved power-flow JSON and the complete parameter workbook. No external
power-flow application or downloaded pickle is needed.

The 50 µs runs take 60,000 steps and the 25 µs runs take 120,000. Each saved
CSV contains 6,001 rows on the same 0.5 ms output grid, including the initial
point. Both disturbances were applied at exactly 1 s. All saved CSV values
are finite. Sampling at 0.5 ms does not resolve every integration step or
establish accuracy of high-frequency switching transients.

The maximum absolute differences between the 50 and 25 µs results, over all
saved samples and all channels in each group, are:

| Quantity | Governor step | Generator trip |
| --- | ---: | ---: |
| Instantaneous phase voltage, pu | 0.0126874 | 4.18823 |
| Three-phase voltage magnitude, pu | 0.000175041 | 4.46649 |
| Rotor speed, pu | 0.0000149218 | 0.00198192 |
| Mechanical power, pu on machine base | 0.000164115 | 0.00995209 |
| Exciter output, pu on exciter base | 0.000196219 | 0.0131627 |
| Stabilizer output, pu | 0.000650735 | 0.0113466 |

The trip's sampled bus-1 voltage magnitude at 1 s grows from **3.595605 pu
to 8.062099 pu** when the time step is halved. This waveform is unsuitable
as a validated transient voltage answer key. Its cause has not been
established; inspect trip reinitialization, network history and numerical
damping before adopting it. The upstream state kernel also continues to
update the disconnected machine's internal states; the trip response plot
therefore shows generators 2 and 3. The CSV preserves generator 1 as well.

The governor step is a better starting point for future controller/network
comparison, but the phase-voltage difference is still about 0.0127 pu, and
two time steps do not prove convergence. These metrics are numerical
differences within ParaEMT, not errors against known truth or GridKit.
[refinement.json](results/refinement.json) also separates the pre-event,
first 20 ms, and later intervals and identifies the maximum-error channel
and time. No pass/fail tolerances have been invented from these results.

Each result directory includes `columns.json` (units), `run.json` (settings,
versions and timing), `run.log`, and `machine_parameters.json` (raw ParaEMT
equivalent-circuit arrays). The finer governor-step run also retains all
saved upstream state vectors in `states.npz`; read using
`numpy.load(..., allow_pickle=False)`. State ordering is documented in
`Initialize.CombineX` in [Lib_BW.py](upstream/Lib_BW.py). Phase voltages use
the peak phase base: multiply by `sqrt(2/3) * V_LL_RMS` to obtain volts.
`vm` is the instantaneous three-phase magnitude, not windowed RMS.

Reproduce from this directory with Python 3.12 (executed with 3.12.3):

```bash
python3 -m venv /tmp/paraemt-reference-venv
/tmp/paraemt-reference-venv/bin/pip install -r upstream/requirements.txt
/tmp/paraemt-reference-venv/bin/python run_reference.py --event governor-step --dt-us 50 --output results/governor_step/dt50us
/tmp/paraemt-reference-venv/bin/python run_reference.py --event governor-step --dt-us 25 --save-states --output results/governor_step/dt25us
/tmp/paraemt-reference-venv/bin/python run_reference.py --event trip --dt-us 50 --output results/trip/dt50us
/tmp/paraemt-reference-venv/bin/python run_reference.py --event trip --dt-us 25 --output results/trip/dt25us
MPLCONFIGDIR=/tmp/paraemt-matplotlib /tmp/paraemt-reference-venv/bin/python plot_results.py
```

The wrapper follows the original driver's step order and history handling;
it selects `systemN=2`, shortens the horizon from 10 to 3 s, exports labeled
compressed CSV instead of pickle, and keeps output spacing fixed while
refining the step. For the governor experiment it moves the already supplied
−0.02 pu reference step from 100 to 1 s and moves the trip beyond the run.
An `ElementTree.getiterator` compatibility alias lets upstream's pinned
`xlrd==1.2.0` read XLSX on current Python. The serial LU path does not need
METIS. [environment.txt](environment.txt) records the installed packages.

Original inputs: `ieee9.raw`, `pfd_9_1_1.json`, and `9bus.xlsx`. The four
short CSV files at this directory's root are worksheet exports, not
simulation results. Original source and license are preserved under
`upstream/` and [LICENSE.md](LICENSE.md). [sources.json](sources.json)
records SHA-256 hashes and upstream Git blob identities. Run
`python3 verify.py` to check the collected artifacts and saved CSV structure;
regenerated floating-point results or timings may have different hashes.

Reference: Xiong et al., *ParaEMT: An Open Source, Parallelizable, and
HPC-Compatible EMT Simulator for Large-Scale IBR-Rich Power Grids*, IEEE
Transactions on Power Delivery 39(2), 911–921 (2024),
[DOI](https://doi.org/10.1109/TPWRD.2023.3342715).
