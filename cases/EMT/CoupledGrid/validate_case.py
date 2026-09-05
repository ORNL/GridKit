#!/usr/bin/env python3
"""Check that the generated EMT case uses the validated coupled line models."""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    case = read_json(HERE / "CoupledGrid.case.json")
    model = read_json(HERE / "lines/model.json")
    fit = read_json(HERE / "lines/fit_validation.json")
    require(hashlib.sha256((HERE / "lines/openline.csv").read_bytes()).hexdigest()
            == fit["source"]["sweep_sha256"], "Raw line data differs from validated sweep")
    devices = {device["id"]: device for device in case["devices"]}
    require(len(devices) == len(case["devices"]), "Duplicate device identifiers")
    buses = {key for key, device in devices.items() if device["class"] == "Bus"}
    lines = [device for device in devices.values() if device["class"] == "LineLumped"]
    require(len(buses) == 14 and len(lines) == 13, "Expected 14 buses and 13 line sections")
    adjacency = {bus: set() for bus in buses}
    for line in lines:
        require(line["submodels"] == model, f'{line["id"]}: stale or modified line coefficients')
        require(line["params"]["N"] == 3 and line["params"]["K"] == 3
                and line["params"]["conductors"] == [1, 2, 3], "Unexpected conductor mapping")
        require(300 <= line["params"]["dx"] <= 600, "Line length is outside validated case range")
        require({"i12a", "i12b", "i12c"}.issubset(line["mon"]), "Missing line-current monitors")
        a, b = (line["inputs"][key] for key in ("bus1", "bus2"))
        require(a in buses and b in buses, "Line endpoint is not a bus")
        adjacency[a].add(b)
        adjacency[b].add(a)
    for device in devices.values():
        if device["class"] == "Switch":
            a, b = (device["inputs"][key] for key in ("bus1", "bus2"))
            adjacency[a].add(b)
            adjacency[b].add(a)
    connected, pending = set(), ["bus_1"]
    while pending:
        bus = pending.pop()
        if bus not in connected:
            connected.add(bus)
            pending.extend(adjacency[bus] - connected)
    require(connected == buses, "A bus remains disconnected when the load-step switch closes")

    # Check the actual serialized realization, including full mutual coupling.
    z = model["Zp"]
    poles = np.array(z["poles"])
    require(len(poles) == 8 and np.all(poles[:, 0] < 0)
            and np.all(poles[:, 1] == 0), "Line poles must be stable and real")
    residues = np.array(z["residues"])
    require(np.all(residues[..., 1] == 0), "Real line poles require real residues")
    foster = residues[..., 0] / poles[:, 0, None, None]
    r0 = np.array(z["D"]) - foster.sum(axis=0)
    for matrix in [r0, np.array(z["E"]), np.array(model["Yp"]["E"]), *foster]:
        require(np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-15), "Line is not reciprocal")
        require(np.linalg.eigvalsh(matrix).min() > -1e-13, "Line realization is not passive")
    require(np.all(np.abs(np.array(z["E"])[~np.eye(3, dtype=bool)]) > 1e-8),
            "Series model is missing mutual inductive coupling")
    require(np.all(np.array(model["Yp"]["E"])[~np.eye(3, dtype=bool)] < 0),
            "Shunt model is missing mutual capacitance")

    require(sum(device["class"] == "Machine" for device in devices.values()) == 2,
            "Expected two synchronous machines")
    for bus in (9, 11, 13):
        pwm, converter, source = (devices[f"{name}_{bus}"] for name in ("pwm", "converter", "filter"))
        require(converter["inputs"]["s"] == pwm["outputs"]["s"], "PWM-converter connection is incomplete")
        require([source["inputs"][name] for name in ("ea", "eb", "ec")]
                == converter["outputs"]["vo"], "Converter-source connection is incomplete")
        require(source["class"] == "DependentVoltageSource" and "Y" in source["submodels"],
                "IBR source lacks its frequency-dependent filter")
    solvers = [read_json(HERE / f"{name}.solver.json") for name in ("low", "middle", "high")]
    expected_mu = [240., np.sqrt(240. * 3600.), 3600.]
    for solver, expected in zip(solvers, expected_mu):
        require(np.isclose(solver["mu"], expected), "Unexpected smoothing parameter")
        require({key: value for key, value in solver.items() if key not in ("mu", "output_file")}
                == {key: value for key, value in solvers[0].items() if key not in ("mu", "output_file")},
                "Mu comparisons use different solver settings or physical case files")
    print(json.dumps({"buses": len(buses), "coupled_line_sections": len(lines),
                      "machines": 2, "PWM_converter_sources": 3,
                      "mu": expected_mu, "status": "passed"}, indent=2))


if __name__ == "__main__":
    main()
