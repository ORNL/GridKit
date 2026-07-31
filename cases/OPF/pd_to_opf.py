#!/usr/bin/env python3
"""Translate a PhasorDynamics case file to an OPF case file.

Network topology, branch impedances, machine ratings, and governor valve
limits are taken from the dynamics case. Fields with no dynamics source are
filled with fictitious values from the constants below.

Generator active power limits come from the governor matched to the machine
through its pmech signal, converted from the governor rating base to the
system base. Machines without a governor get pmin 0 and pmax at rating.

Usage: pd_to_opf.py <input.case.json> <output.opf.json>
"""

import json
import sys

# Fictitious values for data a dynamics case does not carry
VMIN = 0.95         # bus voltage limits [pu]
VMAX = 1.05
QFRAC = 0.5         # generator q limits as a fraction of the p.u. mva rating
COST_C0 = 0.0
COST_C1 = 10.0      # linear cost, stepped per generator so units are not identical
COST_C1_STEP = 2.0
COST_C2 = 0.01
COST_C1_GRID = 25.0  # linear cost of a grid equivalent at an infinite bus

MACHINE_CLASSES = ("Genrou", "Gensal", "GenClassical")
LOAD_CLASSES = ("Load", "LoadZ", "LoadZIP")
GOVERNOR_CLASSES = ("Tgov1",)
TRATE_DEFAULT = 100.0  # MVA, mirrors the Tgov1 component default


def translate(case):
    header = case["header"]
    base_mva = header["va_base"] / 1e6

    machines = [d for d in case["devices"] if d["class"] in MACHINE_CLASSES]
    loads = [d for d in case["devices"] if d["class"] in LOAD_CLASSES]
    branches = [d for d in case["devices"] if d["class"] == "Branch"]

    # Governors are matched to machines through the shared pmech signal
    governors = {}
    for governor in case["devices"]:
        if governor["class"] in GOVERNOR_CLASSES:
            signal = governor["ports"].get("pmech")
            if signal is not None:
                governors[signal] = governor

    # The slack bus is the first infinite bus if present, otherwise the bus
    # of the machine with the largest mva rating.
    infinite = [b["number"] for b in case["buses"] if b["class"] == "infinite_bus"]
    if infinite:
        slack = infinite[0]
    else:
        slack = max(machines, key=lambda m: m["params"].get("mva", base_mva))["ports"]["bus"]

    # Bus names are not always unique in dynamics cases, so colliding names
    # are suffixed with the bus number
    name_counts = {}
    for bus in case["buses"]:
        name_counts[bus["name"]] = name_counts.get(bus["name"], 0) + 1

    buses = []
    for bus in sorted(case["buses"], key=lambda b: b["number"]):
        unique = name_counts[bus["name"]] == 1
        buses.append({
            "class": "Slack" if bus["number"] == slack else "Bus",
            "id": bus["name"] if unique else f'{bus["name"]}_{bus["number"]}',
            "params": {
                "number": bus["number"],
                "kv": bus["params"]["kv"],
                "vmin": VMIN,
                "vmax": VMAX,
            },
        })

    devices = []
    for branch in branches:
        devices.append({
            "class": "Branch",
            "buses": {"from": branch["ports"]["bus1"], "to": branch["ports"]["bus2"]},
            "id": branch["id"],
            "params": {key: branch["params"].get(key, 0.0) for key in ("R", "X", "G", "B")},
        })

    for index, machine in enumerate(machines):
        mva = machine["params"].get("mva", base_mva)
        rating = round(mva / base_mva, 6)

        governor = governors.get(machine["ports"].get("pmech"))
        if governor is None:
            pmin, pmax = 0.0, rating
        else:
            trate = governor["params"].get("Trate", TRATE_DEFAULT)
            pmin = round(governor["params"].get("Pvmin", 0.0) * trate / base_mva, 6)
            pmax = round(governor["params"].get("Pvmax", 1.0) * trate / base_mva, 6)

        devices.append({
            "class": "Generator",
            "buses": {"bus": machine["ports"]["bus"]},
            "id": machine["id"],
            "params": {
                "pmin": pmin,
                "pmax": pmax,
                "qmin": round(-QFRAC * rating, 6),
                "qmax": round(QFRAC * rating, 6),
                "mva": mva,
                "c0": COST_C0,
                "c1": COST_C1 + COST_C1_STEP * (index % 6),
                "c2": COST_C2,
            },
        })

    # Each infinite bus becomes a grid equivalent generator
    for index, number in enumerate(infinite):
        devices.append({
            "class": "Generator",
            "buses": {"bus": number},
            "id": f"EQ{index + 1}",
            "params": {
                "pmin": -10.0,
                "pmax": 10.0,
                "qmin": -10.0,
                "qmax": 10.0,
                "mva": 10.0 * base_mva,
                "c0": 0.0,
                "c1": COST_C1_GRID,
                "c2": 0.0,
            },
        })

    for load in loads:
        devices.append({
            "class": "Load",
            "buses": {"bus": load["ports"]["bus"]},
            "id": load["id"],
            "params": {},
        })

    return {
        "header": {
            "format_version": 0,
            "format_revision": 1,
            "case_name": header["case_name"],
            "case_description": "Derived from the PhasorDynamics case with fictitious data where no dynamics source exists",
        },
        "params": {
            "freq_base": header["freq_base"],
            "va_base": header["va_base"],
        },
        "buses": buses,
        "devices": devices,
    }


def validate(opf):
    numbers = {bus["params"]["number"] for bus in opf["buses"]}
    for name, items in (("bus", opf["buses"]), ("device", opf["devices"])):
        ids = [item["id"] for item in items]
        duplicates = {i for i in ids if ids.count(i) > 1}
        if duplicates:
            raise SystemExit(f"Duplicate {name} ids: {sorted(duplicates)}")
    for device in opf["devices"]:
        for number in device["buses"].values():
            if number not in numbers:
                raise SystemExit(f"Device {device['id']} references unknown bus {number}")
    slacks = [b for b in opf["buses"] if b["class"] == "Slack"]
    if len(slacks) != 1:
        raise SystemExit(f"Expected one slack bus, found {len(slacks)}")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__.strip())

    with open(sys.argv[1]) as stream:
        case = json.load(stream)

    opf = translate(case)
    validate(opf)

    with open(sys.argv[2], "w") as stream:
        json.dump(opf, stream, indent=4)
        stream.write("\n")

    counts = {c: sum(1 for d in opf["devices"] if d["class"] == c)
              for c in ("Branch", "Generator", "Load")}
    print(f"{sys.argv[2]}: {len(opf['buses'])} buses, " +
          ", ".join(f"{n} {c}" for c, n in counts.items()))


if __name__ == "__main__":
    main()
