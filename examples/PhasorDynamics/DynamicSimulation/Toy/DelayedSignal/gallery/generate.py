"""Reproduce the gallery using the configured GridKit Makefiles build."""
import csv
import hashlib
import io
import json
import re
from pathlib import Path
import shlex
import subprocess
import tempfile
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = next(path for path in HERE.parents if (path / "GridKit/Solver").is_dir())
BUILD = REPO / "build"


def read_data(name, kind="samples"):
    with zipfile.ZipFile(HERE / "data.zip") as archive:
        return np.genfromtxt(
            io.BytesIO(archive.read(f"{name}_{kind}.csv")), delimiter=",", names=True
        )


def generate():
    link_dir = BUILD / "tests/IntegrationTests/PhasorDynamics"
    flags = (BUILD / "application/PhasorDynamics/CMakeFiles/DynamicSimulation.dir/flags.make").read_text()
    includes = shlex.split(
        next(line.split("=", 1)[1] for line in flags.splitlines()
             if line.startswith("CXX_INCLUDES"))
    )
    link = shlex.split(
        (link_dir / "CMakeFiles/test_pd_signal_delay.dir/link.txt").read_text()
    )
    libraries = link[link.index("-o") + 2:]
    with tempfile.TemporaryDirectory(prefix="gridkit-delay-gallery-") as directory:
        temporary = Path(directory)
        executable = temporary / "gallery"
        subprocess.run(
            [link[0], "-std=c++20", "-O2", "-fvisibility=hidden", *includes, str(HERE / "gallery.cpp"),
             "-o", str(executable), *libraries], cwd=link_dir, check=True
        )
        with zipfile.ZipFile(
            HERE / "data.zip", "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
        ) as archive:
            for run in json.loads((HERE / "runs.json").read_text()):
                (temporary / "run.json").write_text(json.dumps(run))
                name = run["name"]
                with (temporary / f"{name}.log").open("w") as log:
                    subprocess.run(
                        [str(executable), "run.json", name], cwd=temporary,
                        stdout=log, stderr=log, check=True
                    )
                for kind in ("samples", "steps", "attempts", "midpoints"):
                    filename = f"{name}_{kind}.csv"
                    archive.write(temporary / filename, filename)
                archive.write(temporary / f"{name}.log", f"{name}.log")
                print(name, flush=True)
    provenance = {
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "source_diff_sha256": hashlib.sha256(subprocess.check_output(
            ["git", "diff", "HEAD", "--", "GridKit", "application", "tests", "buildsystem/spack_repo/gridkit/packages/sundials"],
            cwd=REPO
        )).hexdigest(),
        "sundials_step_check_patch_sha256": hashlib.sha256(
            (REPO / "buildsystem/spack_repo/gridkit/packages/sundials/idas-step-check.patch").read_bytes()
        ).hexdigest(),
        "compiler": link[0],
        "description": "Smooth manufactured sources linked to build/ GridKit libraries",
        "error_samples": "Uniform monitor times plus each accepted step midpoint",
        "prehistory": "Analytic input function for t < 0",
    }
    (HERE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def summarize():
    rows = []
    for run in json.loads((HERE / "runs.json").read_text()):
        name = run["name"]
        samples = read_data(name)
        steps = read_data(name, "steps")
        attempts = read_data(name, "attempts")
        midpoints = read_data(name, "midpoints")
        for values in (samples, steps, attempts, midpoints):
            assert all(np.isfinite(values[key]).all() for key in values.dtype.names)
        assert np.all(np.diff(samples["t"]) > 0)
        assert np.allclose(steps["h"], np.diff(np.r_[0, steps["t"]]))
        assert np.allclose(steps["h_over_tau"], steps["h"] / run["tau"])
        assert np.all(steps["error_estimate"] <= 1)
        assert abs(steps["t"][-1] - run["tf"]) < 1e-12
        error = np.r_[abs(samples["z"] - samples["exact_z"]),
                      abs(midpoints["z"] - midpoints["exact_z"])]
        assert max(error) < 20 * run["rel_tol"], (name, float(max(error)))
        with zipfile.ZipFile(HERE / "data.zip") as archive:
            log = archive.read(f"{name}.log").decode()
        native_errors = int(re.search(r"native_error_failures=(\d+)", log)[1])
        nonlinear_failures = int(re.search(r"nonlinear_failures=(\d+)", log)[1])
        rows.append({
            **run,
            "accepted": len(steps),
            "rejected": int(sum(attempts["rejected"])),
            "native_error_failures": native_errors,
            "nonlinear_failures": nonlinear_failures,
            "max_abs_z_error": float(max(error)),
            "max_endpoint_z_error": float(max(abs(steps["z"] - steps["exact_z"]))),
            "min_h_tau": float(min(steps["h_over_tau"])),
            "median_h_tau": float(np.median(steps["h_over_tau"])),
            "max_h_tau": float(max(steps["h_over_tau"])),
            "peak_segments": int(max(steps["segments"])),
            "last_segments": int(steps["segments"][-1]),
        })
    with (HERE / "summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (HERE / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    generate()
    summarize()
    subprocess.run(["python3", "-B", str(HERE / "plot_gallery.py")], check=True)
