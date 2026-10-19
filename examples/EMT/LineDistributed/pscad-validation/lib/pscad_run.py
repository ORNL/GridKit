"""Run PSCAD projects in one PSCAD session and export their recordings as CSV.

Follows GridWorkbench's examples/pscad/study.py: each project is copied to a
work folder (PSCAD builds beside its project), run with PSOUT output, and the
output channels are read back by their stable ``Group`` names with mhi.psout.
"""

import csv
import os
import shutil
import time
from pathlib import Path

import numpy as np

CHANNELS = ["t"] + [f"{q}_{p}" for q in ("e", "v1", "v2", "is", "i2", "rf") for p in "abc"]


class Session:
    """One PSCAD application, launched on first use."""

    def __init__(self, timeout=120):
        self._pscad = None
        self._timeout = timeout

    @property
    def pscad(self):
        if self._pscad is None:
            import mhi.pscad

            os.environ.pop("NoDefaultCurrentDirectoryInExePath", None)
            self._pscad = mhi.pscad.launch(version="5.1.0", silence=True, minimize=True, timeout=self._timeout)
        return self._pscad

    def close(self):
        if self._pscad is not None:
            try:
                self._pscad.quit()
            finally:
                self._pscad = None

    def simulate(self, project: Path, work: Path, results: Path):
        """Run ``project`` from a copy in ``work``; return the PSOUT copied into ``results``."""
        project = Path(project)
        work = Path(work) / project.stem
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        copy = work / project.name
        shutil.copy2(project, copy)
        pscad = self.pscad
        pscad.load(str(copy))
        run = pscad.project(project.stem)
        started = time.time()
        try:
            run.parameters(PlotType="PSOUT")
            run.run()
            for _ in range(40):
                if "Simulation has ended. Status code = 0" in run.output():
                    break
                time.sleep(0.5)
            else:
                problems = [m.text for m in run.messages() if m.status not in ("normal", "warning")]
                raise RuntimeError("\n".join([run.output(), *problems]))
            source = Path(run.temp_folder) / f"{project.stem}.psout"
            if not source.is_file():
                raise FileNotFoundError(f"PSCAD did not produce {source.name}")
            results.mkdir(parents=True, exist_ok=True)
            artifact = results / source.name
            shutil.copy2(source, artifact)
            # Keep the line-constants outputs (fitted data and logs) beside the recording.
            for extra in Path(run.temp_folder).glob("*.tlo"):
                shutil.copy2(extra, results / extra.name)
            for extra in Path(run.temp_folder).glob("*.log"):
                shutil.copy2(extra, results / (project.stem + "_" + extra.name))
            for extra in Path(run.temp_folder).glob("*.out"):
                shutil.copy2(extra, results / extra.name)
            messages = [f"[{m.status}] {m.text}" for m in run.messages()]
            (results / (project.stem + ".messages.txt")).write_text("\n".join(messages), encoding="utf-8")
            return artifact, time.time() - started
        except Exception:
            results.mkdir(parents=True, exist_ok=True)
            (results / (project.stem + ".failure.log")).write_text(run.output(), encoding="utf-8")
            build = Path(run.temp_folder)
            if build.is_dir():
                shutil.copytree(build, results / (project.stem + "_failed_build"), dirs_exist_ok=True)
            raise
        finally:
            try:
                run.unload()
            except Exception:
                pass


def series(artifact: Path) -> dict:
    """Recorded channels of a single-run PSOUT by their Group names: (domain, values)."""
    from mhi.psout import File

    channels = {}
    with File(str(artifact)) as output:
        if output.num_runs != 1:
            raise ValueError(f"{artifact} contains {output.num_runs} runs")
        for trace, path in output.run(0).trace_paths("**/*[@Source='PGB']/Record/*"):
            name = trace.variables()["Group"]
            if name in channels:
                raise ValueError(f"channel {name!r} is ambiguous: {path}")
            channels[name] = (np.array(trace.domain.data, copy=True), np.array(trace.data, copy=True))
    return channels


def export_csv(artifact: Path, destination: Path, channels=CHANNELS):
    """Write the recorded channels on the study clock to ``destination``."""
    data = series(artifact)
    missing = [name for name in channels if name not in data]
    if missing:
        raise KeyError(f"{artifact.name}: missing channels {missing}; have {sorted(data)}")
    domain = data["t"][0]
    for name in channels:
        if not np.array_equal(data[name][0], domain):
            raise ValueError(f"{artifact.name}: {name} is not on the common sample grid")
    columns = np.column_stack([data[name][1] for name in channels])
    np.savetxt(destination, columns, delimiter=",", fmt="%.10e", comments="", header=",".join(channels))
    return destination
