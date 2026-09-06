# Requirements for an equivalent GridKit case

Assessed on `lukel/emt-playground`, GridKit revision
`4bd2812605dd0bbbc1352668bc750497c17347f3`. The obstacle is implementing and
wiring the missing EMT controllers. This collection contains no runnable
GridKit case and makes no GridKit/ParaEMT agreement claim.

| Part | ParaEMT case | Current GridKit EMT status |
| --- | --- | --- |
| Machines | Three GENROU equivalent circuits with two d-axis and two q-axis rotor circuits | `Machine` has corresponding flux/current structure; parameter conversion below is feasible but unvalidated |
| Excitation | SEXS, with terminal-voltage measurement and IEEEST input | EMT exposes IEEET1; SEXS-PTI exists only in PhasorDynamics |
| Governors | GAST on all three machines | EMT exposes TGOV1; GASTPTI exists only in PhasorDynamics |
| Stabilizers | IEEEST on all three machines | No EMT IEEEST; PhasorDynamics implementation exists |
| Network | Six pi lines and three series transformer impedances on per-unit voltage bases | `LineLumped` can represent the continuous pi/RL circuits after referring voltages to a common physical base |
| Loads | Three constant series RL loads for this case and load option 1 | `LoadZ` can represent the continuous RL loads |
| Events | Governor reference step or generator disconnection at 1 s | EMT application schedules switch open/close; no governor-reference-step event is exposed |
| Numerical method | Partitioned updates with explicit control steps, damped companion circuits and trip reinitialization | IDA DAE integration; numerical damping and event limits must be assessed separately from physical parameter equivalence |

The available EMT device data and registration can be checked in
[ContainerData.hpp](../../../GridKit/Model/EMT/ContainerData.hpp) and
[ContainerDataJSONParser.hpp](../../../GridKit/Model/EMT/ContainerDataJSONParser.hpp).
The existing phasor implementations are
[SEXS-PTI](../../../GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/),
[GASTPTI](../../../GridKit/Model/PhasorDynamics/Governor/GASTPTI/) and
[IEEEST](../../../GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/).
They inherit the PhasorDynamics component interface and cannot be named as
EMT devices in the current builder. Porting or adapting them requires EMT
signals/ports, initialization, parser/factory registration and Jacobian
support. Substituting TGOV1/IEEET1, constant mechanical power/field voltage,
or removing the stabilizers would define a different experiment.

The missing controls are active, not unused workbook metadata.
`Initialize.InitExc/InitGov/InitPss` initialize them, `CombineX` allocates
their states, and `EmtSimu.updateX` passes their data to `numba_updateX`.
The saved `gen*_pm_pu`, `gen*_efd_pu` and `gen*_pss_vs_pu` trajectories show
their response. Match the actual pinned kernels, including limiter and
zero-time-constant behavior, when adapting existing GridKit classes.

## Machine parameter mapping

Do not copy GENROU reactances/time constants into similarly named EMT
fundamental parameters without conversion. ParaEMT's
`DyData.ToEquiCirData` and `Initialize.MergeMacG` in
[Lib_BW.py](upstream/Lib_BW.py) expose the actual winding matrix used in
simulation. The saved `machine_parameters.json` contains those raw `ec_*`
arrays. GridKit's machine equations are in
[MachineImpl.hpp](../../../GridKit/Model/EMT/Component/Source/Machine/MachineImpl.hpp).

Comparing the continuous flux equations gives this candidate mapping:

| GridKit parameter | ParaEMT expression |
| --- | --- |
| `Ll`, `Lmd`, `Lmq`, `L0`, `Rs` | `ec_Ll`, `ec_Lad`, `ec_Laq`, `ec_L0`, `ec_Ra` |
| `Llfd` | `ec_Lffd - ec_Lad` |
| `Ll1d` | `ec_L11d - ec_Lad` |
| `Ll1q` | `ec_L11q - ec_Laq` |
| `Ll2q` | `ec_L22q - ec_Laq` |
| `Rfd`, `R1d`, `R1q`, `R2q` | Corresponding `ec_R* / ws` |
| `S`, `V`, `f` | Machine MVA base × 1e6, nominal line-line RMS kV × 1e3, 60 Hz |
| `p0`, `q0` | Solved generator MW/Mvar × 1e6 |
| `H`, `F`, `S10`, `S12` | Workbook `H`, 0, 0, 0 for this case |

The resistance scaling follows GridKit's rotor equation
`(1/ws) * d(psi_fd)/dt + Rfd * ifd - efd = 0`; ParaEMT's rotor circuit
uses an unscaled flux derivative. Its exciter conversion changes by the
same factor. The d/q axis angle and current conventions still need to be
aligned before comparing internal states.

There is a subtle difference between intermediate and assembled ParaEMT
parameters: `ec_L1d = 0.2`, but the assembled d-axis damper self-inductance
is `ec_L11d = 1.4083333333`. With `ec_Lad = 1.3`, GridKit's candidate
`Ll1d` is **0.1083333333**, not 0.2. Similarly, `Ll2q` should be
`1.3586956522 - 1.25 = 0.1086956522`, not intermediate `ec_L2q = 0.125`.
These expressions are an equation-based mapping proposal, not a validated
GridKit machine initialization or completed case conversion.

## Network and measurements

The pinned ParaEMT network kernel in [lib_numba.py](upstream/lib_numba.py)
uses dimensionless bus voltages and series RL transformer branches; it
does not introduce explicit ideal winding-ratio equations. The supplied
transformer taps are all one on their respective voltage bases. An
equivalent continuous network can therefore be referred to, for example,
230 kV everywhere, with 100 MVA network base and `Z_base = 529 ohm`.
Machine voltage ratings and bus initial voltages must be referred together;
keep each machine's own MVA base and per-unit internal parameters.

On that base, use diagonal phase matrices with line
`R = Re(Z_pu) * Z_base`, `L = Im(Z_pu) * Z_base / ws` and total pi
`C = B_pu / (ws * Z_base)` (half at each end). Transformer branches use
their series impedance with zero shunt. Loads use
`Z_pu = Vm^2 / conjugate((P+jQ)/100 MVA)` and the same R/L conversion.
This is a representation of the case's balanced continuous network; it
does not establish transformer zero-sequence, saturation or unbalanced
fault equivalence to a physical winding model.

ParaEMT also adds time-step-dependent damping to line inductors/capacitors
and inductive loads (`damptrap`, `Rp`, `Rs` in `numba_InitNet`) and uses
`Init_mac_alpha = 99/101`. Copying the physical RLC parameters alone will
not reproduce finite-step damping exactly. First compare initialization,
60 Hz phase relationships and the governor step; then refine the time step
and GridKit tolerance separately. The supplied trip spike is a reason to
investigate event behavior before setting waveform acceptance tolerances.

The SEXS input uses the filtered voltage magnitude state from ParaEMT's
bus measurement block; its time constant is in the workbook's `vm` sheet.
IEEEST uses rotor speed deviation. Preserve those signal definitions,
their initialization, the 150/250/100 MVA machine bases and generator
terminal voltage bases of 16.5/18/13.8 kV. The solved JSON, rather than a
generic IEEE 9-bus operating point, is the reference initialization.

After these model and event gaps are addressed, author the EMT `devices`,
`inputs`, `outputs`, state and solver JSON, then compare physical bus
waveforms, rotor speed, mechanical power, exciter output and stabilizer
output against the saved columns. A simulation and tolerance study remain
necessary before calling that future case equivalent or validated.
