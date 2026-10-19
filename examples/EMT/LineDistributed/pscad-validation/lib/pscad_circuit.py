"""PSCAD validation projects built with GridWorkbench's PSCAD authoring internals.

Each project is the LineDistributed demonstration circuit in PSCAD:

    three-phase controlled source e(t) -- Rs (+Ls) -- bus 1 == line == bus 2 -- R(-L) load
                                                                           |
                                                              variable fault resistance

A scripted ``Forcing`` component evaluates the same source and fault functions
as the GridKit example (finite C2 transitions, additive pulses, amplitude
events), so both tools integrate identical forcing. The line is either a
frequency-dependent phase model from the tower geometry or a Bergeron line
from 60 Hz sequence data. Instantaneous terminal voltages and currents are
recorded in SI units on PSCAD's study clock (EMTDC time less the start-up).
"""

from dataclasses import dataclass
from typing import ClassVar, Mapping, Optional

import numpy as np
from gridworkbench.io.pscad import parts as p
from gridworkbench.io.pscad.circuit import Circuit
from gridworkbench.io.pscad.compile import build
from gridworkbench.io.pscad.core import PHASE, THREE_PHASE, Net, Part, Port, monitor, name_ref, omit_none, param
from gridworkbench.io.pscad.native import encode
from gridworkbench.io.pscad.project import Blocks, LineResource, Page, Project, Run

KILO = 1e3
NEGLIGIBLE_R = 1e-12  # pu/m: PSCAD's manual-data Bergeron line with no series loss


# --- master-library parts the GridWorkbench adapter does not declare -------------------


@dataclass(frozen=True)
class VariableResistor(Part):
    """Variable resistance (ohm) read from a signal every time step."""

    DEFINITION: ClassVar[str] = "varrlc"
    PORTS: ClassVar[Mapping[str, Port]] = {
        "A": Port("electrical", (0, 0)),
        "B": Port("electrical", (2, 0)),
    }
    name: Optional[str] = name_ref()
    branch: int = param("RLC", 0, init=False)
    resistance: Optional[str] = omit_none("R", signal="input")
    rate_effects: int = param("dLdC", 0, init=False)


@dataclass(frozen=True)
class Ammeter(Part):
    """Branch current in kA, positive from N1 to N2."""

    DEFINITION: ClassVar[str] = "ammeter"
    PORTS: ClassVar[Mapping[str, Port]] = {
        "N1": Port("electrical", (0, 0)),
        "N2": Port("electrical", (2, 0)),
    }
    name: Optional[str] = omit_none("Name", signal="output")


def forcing_part(code):
    """A scripted component class whose Fortran body is ``code``: study time in, the
    three source voltages (kV) and the three fault resistances (ohm) out."""

    @dataclass(frozen=True)
    class Forcing(Part):
        LIBRARY: ClassVar[str] = ""
        DEFINITION: ClassVar[str] = "Forcing"
        clock: Optional[str] = omit_none("Clock", signal="input")
        Ea: Optional[str] = monitor("Ea")
        Eb: Optional[str] = monitor("Eb")
        Ec: Optional[str] = monitor("Ec")
        Ra: Optional[str] = monitor("Ra")
        Rb: Optional[str] = monitor("Rb")
        Rc: Optional[str] = monitor("Rc")
        SCRIPTS: ClassVar[tuple] = (("Fortran", code),)

    return Forcing


# --- line definitions --------------------------------------------------------------------


@dataclass(frozen=True)
class TowerLine(p.GeometryLine):
    """A frequency-dependent phase-model line from a universal tower with bundles and
    eliminated ground wires, and explicit curve-fitting options."""

    sub_conductors: int = 1
    bundle_spacing: float = 0.4572
    ground_wires: tuple = ()  # (x, h, radius, dc_resistance_ohm_per_km) per wire
    shunt_conductance: float = 1e-11
    fit_start_hz: float = 0.01
    fit_end_hz: float = 1e6
    fit_points: int = 1000
    yc_poles: int = 20
    yc_error_percent: float = 0.05
    h_poles: int = 20
    h_error_percent: float = 0.05
    dc_correction: bool = True
    passivity_scan: bool = True
    detailed_output: bool = True


_original_write_line = encode._write_line


def _write_line(project, line):
    if not isinstance(line, TowerLine):
        return _original_write_line(project, line)
    definition = project.definition.create_tline(line.name)
    options = definition.freq_phase_options.parameters
    options.Interp1 = True
    options.Output = bool(line.detailed_output)
    options.FS = f"{line.fit_start_hz} [Hz]"
    options.FE = f"{line.fit_end_hz} [Hz]"
    options.Numf = int(line.fit_points)
    options.YMaxP = line.yc_poles
    options.YMaxE = f"{line.yc_error_percent} [%]"
    options.AMaxP = line.h_poles
    options.AMaxE = f"{line.h_error_percent} [%]"
    options.DCenab = bool(line.dc_correction)
    options.CPASS = "Scan only" if line.passivity_scan else "Disable"
    definition.ground.parameters.GRRES = f"{line.earth_resistivity} [ohm*m]"
    radius, inner = line.conductors[0][2:]
    ground = {}
    if line.ground_wires:
        x, h, g_radius, g_resistance = line.ground_wires[0]
        ground = dict(
            GDataType="direct",
            GWSame=True,
            ElimGW=True,
            RadiusG=f"{g_radius} [m]",
            DCResG=f"{g_resistance} [ohm/km]",
            PERMG=1.0,
            SAGG="0.0 [m]",
        )
    definition.tower_universal(
        0.0,
        [(x, y) for x, y, _, _ in line.conductors],
        [(x, h) for x, h, _, _ in line.ground_wires],
        TName=line.name,
        Transp=True,
        STL="No, all together",
        ShuntG=f"{line.shunt_conductance} [mho/m]",
        SHSag=False,
        CDataType="direct",
        Hollow="hollow" if inner > 0 else "solid core",
        RadiusC=f"{radius} [m]",
        RadiusCin=f"{inner} [m]",
        DCResC=f"{line.dc_resistance} [ohm/km]",
        PERMC=1.0,
        SAGC="0.0 [m]",
        NCondB=line.sub_conductors,
        BundSym="symmetrical",
        BSP=f"{line.bundle_spacing} [m]",
        **ground,
    )
    return definition


encode._write_line = _write_line


def tower_line(name, tower, frequency_hz=60.0, **options):
    """A :class:`TowerLine` from a ``tower.json`` description."""
    conductor = tower["conductor"]
    return TowerLine(
        name,
        tower["length_m"],
        frequency_hz,
        tuple((c["x"], c["h"], conductor["radius_m"], conductor["inner_radius_m"]) for c in tower["conductors"]),
        dc_resistance=conductor["dc_resistance_ohm_per_km"],
        earth_resistivity=tower["earth_resistivity_ohm_m"],
        sub_conductors=tower["sub_conductors"],
        bundle_spacing=tower["bundle_spacing_m"] or 0.4572,
        ground_wires=tuple((g["x"], g["h"], g["radius"], g["dc_resistance_ohm_per_km"]) for g in tower["ground_wires"]),
        shunt_conductance=tower["shunt_conductance_S_per_m"],
        **options,
    )


def bergeron_line(name, sequence, frequency_hz=60.0, mva=100.0):
    """A lossless PSCAD Bergeron line from the 60 Hz sequence L and C."""
    kv = sequence["voltage_kv"]
    zb = kv**2 / mva
    omega = 2 * np.pi * frequency_hz
    return p.BergeronLine(
        name,
        sequence["length_m"],
        frequency_hz,
        kv,
        mva,
        r1=NEGLIGIBLE_R,
        x1=omega * sequence["L1"] / zb,
        b1=omega * sequence["C1"] * zb,
        r0=NEGLIGIBLE_R,
        x0=omega * sequence["L0"] / zb,
        b0=omega * sequence["C0"] * zb,
    )


# --- forcing script ----------------------------------------------------------------------


def fortran(value):
    """A double-precision Fortran literal."""
    return f"{float(value):.16E}".replace("E", "D")


RAMP = "X*X*X*(1.0D1 + X*(-1.5D1 + 6.0D0*X))"


def forcing_code(spec, startup):
    """Fixed-form Fortran evaluating the GridKit example's source and fault functions.

    ``spec`` is one GridKit run (source, initialization, events); ``startup`` is the
    EMTDC start-up before study time zero, during which a sinusoidal initial state is
    ramped in smoothly over 50 ms.
    """
    lines = ["#LOCAL REAL TS", "#LOCAL REAL X", "#LOCAL REAL ENV", "#LOCAL REAL W"]
    for name in ("AMP", "PLS", "G", "VPK", "PH"):
        lines += [f"#LOCAL REAL {name}{phase}" for phase in "ABC"]
    source = spec["source"]
    steady = spec["initialization"] == "sinusoidal"
    body = ["      TS = $Clock", f"      W = {fortran(2 * np.pi * source['frequency'])}"]
    for k, phase in enumerate("ABC"):
        body += [
            f"      VPK{phase} = {fortran(np.sqrt(2) * source['rms'][k])}",
            f"      PH{phase} = {fortran(source['phase'][k])}",
            f"      PLS{phase} = 0.0D0",
            f"      G{phase} = 0.0D0",
        ]
    if steady and startup > 0:
        body += [
            f"      X = (TS + {fortran(startup)})/5.0D-2",
            "      IF (X .LT. 0.0D0) X = 0.0D0",
            "      IF (X .GT. 1.0D0) X = 1.0D0",
            f"      ENV = {RAMP}",
        ] + [f"      AMP{phase} = ENV" for phase in "ABC"]
    else:
        body += [f"      AMP{phase} = {fortran(1.0 if steady else 0.0)}" for phase in "ABC"]
    target = {"source": "AMP", "voltage": "PLS", "fault": "G"}
    for event in spec["events"]:
        body += [
            f"      X = (TS - {fortran(event['on'])})/{fortran(event['rise'])}",
            "      IF (X .LT. 0.0D0) X = 0.0D0",
            "      IF (X .GT. 1.0D0) X = 1.0D0",
            f"      ENV = {RAMP}",
            f"      X = (TS - {fortran(event['off'])})/{fortran(event['fall'])}",
            "      IF (X .LT. 0.0D0) X = 0.0D0",
            "      IF (X .GT. 1.0D0) X = 1.0D0",
            f"      ENV = ENV - {RAMP}",
        ]
        for phase in event["phases"]:
            name = target[event["type"]] + "ABC"[phase - 1]
            body.append(f"      {name} = {name} + {fortran(event['value'])}*ENV")
    for phase in "ABC":
        body += [
            f"      $E{phase.lower()} = PLS{phase} + AMP{phase}*VPK{phase}*COS(W*TS + PH{phase})",
            f"      $E{phase.lower()} = $E{phase.lower()}*1.0D-3",
            f"      $R{phase.lower()} = 1.0D0/MAX(G{phase}, 1.0D-9)",
        ]
    for line in body:
        if len(line) > 72:
            raise ValueError(f"Fortran line exceeds 72 columns: {line}")
    return "\n".join(lines + body) + "\n"


# --- the circuit ---------------------------------------------------------------------------


def project(name, title, spec, line_definition, *, startup, step, sample, description=""):
    """The validation project for one GridKit run ``spec`` on ``line_definition``."""
    c = Circuit(("Main",))
    source, load = spec["source"], spec["load"]
    clock = c.net("clock")
    c.add("clock", p.SimulationTime(), OUT=clock)
    study_time = c.sub(clock, float(startup), name="study_time")
    c.probe(study_time, name="t", group="t", unit="s")
    forcing = {key: c.net(key) for key in ("Ea", "Eb", "Ec", "Ra", "Rb", "Rc")}
    c.add("forcing", forcing_part(forcing_code(spec, startup))(), clock=study_time, **forcing)
    ground = c.net("ground", PHASE)
    c.add("ground", p.Ground(), A=ground)
    bus1 = {phase: c.net(f"bus1_{phase}", PHASE) for phase in "abc"}
    bus2 = {phase: c.net(f"bus2_{phase}", PHASE) for phase in "abc"}
    line_end = {phase: c.net(f"line2_{phase}", PHASE) for phase in "abc"}
    for k, phase in enumerate("abc"):
        e = forcing["E" + phase]
        c.probe(e, name=f"e_{phase}", group=f"e_{phase}", scale=KILO, unit="V")
        terminal = c.net(f"source_{phase}", PHASE)
        # source_1 imposes V(NA) - V(NB) = Mag: the network side is NA, the return NB.
        c.add(f"source_{phase}", p.ControlledVoltageSource(name=f"E{phase}"), NA=terminal, NB=ground, Mag=e)
        current = c.net(f"is_{phase}")
        inductance = source.get("inductance", [0, 0, 0])[k]
        c.add(
            f"impedance_{phase}",
            p.SeriesRlcBranch(
                name=f"Zs{phase}",
                resistance_enabled=1,
                inductance_enabled=int(inductance > 0),
                capacitance_enabled=0,
                resistance=source["resistance"][k],
                inductance=inductance if inductance > 0 else None,
                current=current,
            ),
            NA=terminal,
            NB=bus1[phase],
        )
        c.probe(current, name=f"is_{phase}", group=f"is_{phase}", scale=KILO, unit="A")
        for end, bus in ((1, bus1), (2, bus2)):
            voltage = c.net(f"v{end}_{phase}")
            c.add(f"meter_v{end}_{phase}", p.NodeLabel(name=f"N{end}{phase}", measure_voltage=1, normalize=0, voltage_signal=voltage), A=bus[phase])
            c.probe(voltage, name=f"v{end}_{phase}", group=f"v{end}_{phase}", scale=KILO, unit="V")
        current2 = c.net(f"i2_{phase}")
        c.add(f"ammeter_{phase}", Ammeter(), N1=line_end[phase], N2=bus2[phase], name=current2)
        c.probe(current2, name=f"i2_{phase}", group=f"i2_{phase}", scale=KILO, unit="A")
        c.add(f"fault_{phase}", VariableResistor(name=f"Rf{phase}"), A=bus2[phase], B=ground, resistance=forcing["R" + phase])
        c.probe(forcing["R" + phase], name=f"rf_{phase}", group=f"rf_{phase}", unit="ohm")
    bus1_abc, bus2_abc, line2_abc = (c.net(key, THREE_PHASE) for key in ("bus1", "bus2", "line2"))
    c.add("bus1_phases", p.Breakout(), N=bus1_abc, N1=bus1["a"], N2=bus1["b"], N3=bus1["c"])
    c.add("line2_phases", p.Breakout(), N=line2_abc, N1=line_end["a"], N2=line_end["b"], N3=line_end["c"])
    c.add("bus2_phases", p.Breakout(), N=bus2_abc, N1=bus2["a"], N2=bus2["b"], N3=bus2["c"])
    c.add("line_send", p.LineEnd(line=line_definition.name, end=1), **{"Z:NC": bus1_abc})
    c.add("line_receive", p.LineEnd(line=line_definition.name, end=2), **{"Z:NC": line2_abc})
    c.add("line", p.LineWire(definition=line_definition.name, length_km=line_definition.length_m / KILO,
                             frequency_hz=line_definition.frequency_hz, remote=True))
    inductive = any(value > 0 for value in load["inductance"])
    c.add(
        "load",
        p.WyeLoad(
            name="Load",
            branch_a=4 if inductive else 1,
            branch_b=4 if inductive else 1,
            branch_c=4 if inductive else 1,
            **{f"resistance_{phase}": load["resistance"][k] for k, phase in enumerate("abc")},
            **({f"inductance_{phase}": load["inductance"][k] for k, phase in enumerate("abc")} if inductive else {}),
        ),
        ABC=bus2_abc,
    )
    page = Page("Main", title, c.freeze(), Blocks())
    return Project(
        name,
        page,
        Run(startup + spec["solver"]["tf"], step, sample, startup),
        (LineResource(("line", line_definition.name), line_definition),),
        description,
    )


def write(project_, path):
    """Compile and write ``project_`` to ``path`` (.pscx)."""
    return encode.write(build(project_), path)
