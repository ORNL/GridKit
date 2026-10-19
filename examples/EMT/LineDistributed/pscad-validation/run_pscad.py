"""Simulate every project in scenarios/pscad with PSCAD and export CSV recordings.

Results go to ``results/pscad/<project>/``: the PSOUT, the CSV on the study clock,
PSCAD's build messages and the line-constants files. Finished projects (a CSV
exists) are skipped unless ``--force`` is given. Requires Windows, PSCAD 5.1 and
the GridWorkbench environment (mhi-pscad, mhi-psout).

    python run_pscad.py [--force] [pattern ...]
"""

import json
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

from pscad_run import Session, export_csv  # noqa: E402

PROJECTS = HERE / "scenarios" / "pscad"
RESULTS = HERE / "results" / "pscad"
WORK = HERE / "results" / "pscad-build"


def main(argv):
    force = "--force" in argv
    patterns = [a for a in argv if not a.startswith("--")]
    projects = sorted(PROJECTS.glob("*.pscx"))
    if patterns:
        projects = [p for p in projects if any(pattern in p.stem for pattern in patterns)]
    session = Session()
    timings = {}
    log = RESULTS / "runs.json"
    if log.exists():
        timings = json.loads(log.read_text())
    failures = []
    try:
        for project in projects:
            results = RESULTS / project.stem
            csv = results / (project.stem + ".csv")
            if csv.exists() and not force:
                print(f"skip {project.stem}", flush=True)
                continue
            started = time.time()
            try:
                artifact, seconds = session.simulate(project, WORK, results)
                export_csv(artifact, csv)
                timings[project.stem] = {"pscad_seconds_including_build": seconds, "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
                print(f"done {project.stem} in {time.time() - started:.1f} s", flush=True)
            except Exception as error:  # keep going; the summary lists failures
                failures.append(project.stem)
                print(f"FAILED {project.stem}: {error}", flush=True)
                traceback.print_exc()
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(json.dumps(timings, indent=2) + "\n")
    finally:
        session.close()
    if failures:
        print("failures:", failures)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
