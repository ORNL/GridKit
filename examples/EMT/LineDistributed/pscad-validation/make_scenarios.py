"""Write the GridKit run suites and the PSCAD projects for every line, model and scenario.

- ``scenarios/scenarios.json``: every expanded run object, for the record.
- ``scenarios/gridkit/<line>__<model>.json``: a ``line_distributed`` suite (defaults and runs).
- ``scenarios/pscad/<project>.pscx``: one PSCAD project per run, plus half-step
  refinements of the frequency-dependent energize, fault and pulse runs.

Run with the GridWorkbench environment after ``make_lines.py``.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

from pscad_circuit import bergeron_line, project, tower_line, write  # noqa: E402
from scenarios import CODES, MODELS, PSCAD_STEPS, STARTUP, project_name, scenarios  # noqa: E402

LINES = HERE / "lines"
GRIDKIT = HERE / "scenarios" / "gridkit"
PSCAD = HERE / "scenarios" / "pscad"
REFINED = {("fd", "energize"), ("fd", "fault"), ("fd", "pulse")}


def main(argv):
    gridkit_only = "--gridkit-only" in argv
    names = [a for a in argv if not a.startswith("--")]
    GRIDKIT.mkdir(parents=True, exist_ok=True)
    PSCAD.mkdir(parents=True, exist_ok=True)
    record = {"models": MODELS, "pscad_steps_s": PSCAD_STEPS, "pscad_startup_s": STARTUP, "lines": {}}
    for line in CODES:
        if names and line not in names:
            continue
        sequence = json.loads((LINES / line / "sequence.json").read_text())
        tower = json.loads((LINES / line / "tower.json").read_text())
        record["lines"][line] = {"sequence": sequence, "runs": {}}
        for model in MODELS:
            runs = scenarios(line, sequence, model)
            record["lines"][line]["runs"][model] = runs
            fit = "fd.fit.json" if model == "fd" else "bergeron.fit.json"
            # The executable's parser accepts only its own fields; the description stays in scenarios.json.
            suite = {"defaults": {"line": f"../../lines/{line}/{fit}"},
                     "runs": [{key: value for key, value in run.items() if key != "description"} for run in runs.values()]}
            (GRIDKIT / f"{line}__{model}.json").write_text(json.dumps(suite, indent=2) + "\n")
            if gridkit_only:
                continue
            definition_name = f"{CODES[line]}_{model}"
            definition = tower_line(definition_name, tower) if model == "fd" else bergeron_line(definition_name, sequence)
            for scenario, spec in runs.items():
                step, sample = PSCAD_STEPS[scenario]
                startup = STARTUP if spec["initialization"] == "sinusoidal" else 0.0
                variants = [(False, step, sample)]
                if (model, scenario) in REFINED:
                    variants.append((True, step / 2, sample))
                for refined, h, s in variants:
                    name = project_name(line, model, scenario, refined)
                    title = f"{sequence['voltage_kv']:g} kV {line}: {MODELS[model].split(' (')[0]} line, {scenario}"
                    built = project(name, title, spec, definition, startup=startup, step=h, sample=s,
                                    description=spec["description"])
                    print(write(built, PSCAD / f"{name}.pscx"), flush=True)
    (HERE / "scenarios" / "scenarios.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
