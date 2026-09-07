"""Run the normal EMT application and compare seeded currents with an RL oracle."""
import argparse
import cmath
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent


def validate(executable, scenario):
    solver_path = ROOT / f"{scenario}.solver.json"
    solver = json.loads(solver_path.read_text())
    case = json.loads((ROOT / solver["system_model_file"]).read_text())
    devices = {device["id"]: device for device in case["devices"][0]["devices"]}
    source = devices["source"]["params"]
    omega = source["omega"]
    seeds = json.loads((ROOT / solver["state_file"]).read_text())["devices"]
    # Absolute bounds for 10 V peak / ampere-scale currents, using 1e-10 IDA
    # tolerances. These cover global integration and interpolation errors.
    current_tolerance = 2e-7
    derivative_tolerance = 2e-5
    maximum = 0.0
    with tempfile.TemporaryDirectory(prefix="gridkit-emt-current-") as directory:
        run = subprocess.run([str(executable), str(solver_path)], cwd=directory,
                             text=True, capture_output=True, check=False)
        if run.returncode:
            raise RuntimeError(run.stdout + run.stderr)
        output = Path(directory)
        layout = json.loads((output / "state.csv.json").read_text())["variables"]
        metadata = {(item["component"], item["local_index"]): item for item in layout}
        with (output / "state.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        expected_rows = round(solver["tmax"] / solver["dt_monitor"]) + 1
        assert len(rows) == expected_rows, (scenario, len(rows), expected_rows)
        for index, row in enumerate(rows):
            time = float(row["time"])
            assert abs(time - index * solver["dt_monitor"]) < 1e-12
            assert all(math.isfinite(float(value)) for value in row.values())
            for name, sign, fields in (("line", 1, "i12"), ("load", -1, "i"), ("resistor", -1, "i")):
                device = devices[name]["params"]
                resistance = device["Rp" if name == "line" else "R"]
                inductance = device.get("Lp" if name == "line" else "L", [[0.0] * 3 for _ in range(3)])
                path = "network." + name
                for phase in range(3):
                    r, induct = resistance[phase][phase], inductance[phase][phase]
                    voltage = math.sqrt(2) * source["E"][phase] * cmath.exp(1j * source["phi"][phase])
                    steady = sign * voltage / complex(r, omega * induct)
                    rotating = steady * cmath.exp(1j * omega * time)
                    initial = seeds.get(path, {}).get(fields + "abc"[phase], 0.0)
                    natural = (initial - steady.real) * math.exp(-r * time / induct) if induct else 0.0
                    expected = rotating.real + natural
                    column = f"{path}[{phase}]"
                    actual = float(row["y:" + column])
                    error = abs(actual - expected)
                    maximum = max(maximum, error)
                    assert error <= current_tolerance, (scenario, time, column, actual, expected)
                    assert metadata[(path, phase)]["differential"] == bool(induct)
                    if induct:
                        derivative = (1j * omega * rotating).real - r / induct * natural
                        assert abs(float(row["yp:" + column]) - derivative) <= derivative_tolerance, (scenario, time, column, "derivative")
                        if index == 0:
                            assert abs(actual - initial) <= 1e-12, (scenario, column, "initial current changed")
        with (output / "mon.csv").open(newline="") as stream:
            monitor = list(csv.reader(stream))
        assert len(monitor) == expected_rows + 1
    print(f"{scenario}: {expected_rows} samples, maximum current error {maximum:.3e} A")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    arguments = parser.parse_args()
    for scenario in ("Steady", "Perturbed", "Default"):
        validate(arguments.exe.resolve(), scenario)
