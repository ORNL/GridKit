"""Summarize PSCAD's frequency-dependent fits from its line-constants logs.

Reads ``results/pscad/V_<code>_fd_fault/V_<code>_fd_fault_<code>_fd.log`` and writes
``lines/<line>/pscad-fit.json`` with the Yc fit (poles, error) and the final
propagation-function attempt (delay groups, poles, delays, errors), beside the
GridKit fit report for the same geometry.
"""

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

from scenarios import CODES, project_name  # noqa: E402

LINES = HERE / "lines"
PSCAD = HERE / "results" / "pscad"


def parse(text):
    yc = re.search(r"Yc:.*?-{8,}\s*\n\s*([\d.]+)\s*%\s+(\d+)", text, re.S)
    attempts = re.split(r"Attempt #\s*\d+", text)[1:]
    groups = []
    phase = None
    if attempts:
        last = attempts[-1]
        for match in re.finditer(r"^\s+(\d+)\s+([\d.]+)\s*%\s+(\d+)\s+([\d.]+)\s*ms\s*$", last, re.M):
            groups.append({"group": int(match.group(1)), "max_error_percent": float(match.group(2)),
                           "poles": int(match.group(3)), "delay_s": float(match.group(4)) * 1e-3})
        hp = re.search(r"Hphase:.*?-{8,}\s*\n\s*([\d.]+)\s*%\s+([\d.]+)\s*%\s+([\d.]+)", last, re.S)
        if hp:
            phase = {"max_error_percent": float(hp.group(1)), "rms_error_percent": float(hp.group(2)), "residue_pole_ratio": float(hp.group(3))}
    return {"Yc": {"max_error_percent": float(yc.group(1)), "poles": int(yc.group(2))} if yc else None,
            "H_attempts": len(attempts), "H_delay_groups": groups, "H_phase_domain": phase}


def main():
    for line, code in CODES.items():
        name = project_name(line, "fd", "fault")
        log = PSCAD / name / f"{name}_{code}_fd.log"
        if not log.exists():
            print("missing", log)
            continue
        summary = parse(log.read_text(errors="replace"))
        gridkit = json.loads((LINES / line / "fit-report.json").read_text())
        summary["gridkit"] = {"Yc_poles": gridkit["Yc"]["poles"], "Yc_max_relative_error": gridkit["Yc"]["max_relative_error"],
                              "H_poles": gridkit["H"]["poles"], "H_max_scaled_error": gridkit["H"]["max_scaled_error"],
                              "delays_s": gridkit["wavefront_delays_s"]}
        (LINES / line / "pscad-fit.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(line, json.dumps(summary))


if __name__ == "__main__":
    main()
