"""Build the gallery driver, record simulations, and redraw diagnostic figures."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = next(path for path in HERE.parents if (path / "GridKit/Solver").is_dir())


def read_data(name, kind="samples"):
    with zipfile.ZipFile(HERE / "data.zip") as archive:
        return np.atleast_1d(np.genfromtxt(
            io.BytesIO(archive.read(f"{name}_{kind}.csv")), delimiter=",", names=True
        ))


def generate(build):
    subprocess.run(
        ["cmake", "--build", str(build), "--target", "signal_delay_gallery", "-j", "10"],
        check=True,
    )
    executable = build / HERE.parent.relative_to(REPO) / "signal_delay_gallery"
    with tempfile.TemporaryDirectory(prefix="gridkit-delay-gallery-") as directory:
        temporary = Path(directory)
        # Publish the complete archive only after every simulation succeeds.
        with zipfile.ZipFile(temporary / "data.zip", "w", zipfile.ZIP_DEFLATED) as archive:
            for run in json.loads((HERE / "runs.json").read_text()):
                (temporary / "run.json").write_text(json.dumps(run))
                name = run["name"]
                with (temporary / f"{name}.log").open("w") as log:
                    subprocess.run([str(executable), "run.json", name], cwd=temporary,
                                   stdout=log, stderr=log, check=True)
                for kind in ("samples", "steps", "attempts", "midpoints", "segments"):
                    filename = f"{name}_{kind}.csv"
                    archive.write(temporary / filename, filename)
                archive.write(temporary / f"{name}.log", f"{name}.log")
                print(name, flush=True)
        (HERE / "data.zip").write_bytes((temporary / "data.zip").read_bytes())
    provenance = {
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "source_diff_sha256": hashlib.sha256(subprocess.check_output(
            ["git", "diff", "HEAD", "--", "GridKit", "application", "tests", "buildsystem",
             str(HERE.parent / "CMakeLists.txt"), str(HERE / "gallery.cpp"),
             str(HERE / "generate.py"), str(HERE / "plot_gallery.py"), str(HERE / "runs.json")],
            cwd=REPO,
        )).hexdigest(),
        "source_files_sha256": {
            str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [REPO / "GridKit/Model/PhasorDynamics/SignalHistoryPolynomial.hpp", HERE / "gallery.cpp"]
        },
        "sundials_step_check_patch_sha256": hashlib.sha256(
            (REPO / "buildsystem/spack_repo/gridkit/packages/sundials/idas-step-check.patch").read_bytes()
        ).hexdigest(),
        "build": str(build),
        "description": "CMake signal_delay_gallery target using the actual GridKit history and IDAS adapter",
        "error_samples": "Uniform monitor times plus each accepted step midpoint",
        "prehistory": "Analytic value and derivative functions",
        "retained_history": "Uncompressed and compressed SignalHistory instances on the same accepted source knots; neither feeds back into the solve",
    }
    (HERE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def summarize():
    rows = []
    for run in json.loads((HERE / "runs.json").read_text()):
        name = run["name"]
        samples, steps, attempts, midpoints = [
            read_data(name, kind) for kind in ("samples", "steps", "attempts", "midpoints")
        ]
        for values in (samples, steps, attempts, midpoints):
            assert all(np.isfinite(values[key]).all() for key in values.dtype.names), name
        assert np.all(np.diff(samples["t"]) > 0), name
        assert np.allclose(steps["h"], np.diff(np.r_[0, steps["t"]])), name
        assert np.allclose(steps["h_over_tau"], steps["h"] / run["tau"]), name
        assert np.all(steps["error_estimate"] <= 1), name
        assert abs(steps["t"][-1] - run["tf"]) < 1e-12, name
        # Independent count: retain each accepted interval whose end >= t - tau.
        retained = np.arange(1, len(steps) + 1) - np.searchsorted(
            steps["t"], steps["t"] - run["tau"], side="left"
        )
        assert np.array_equal(steps["uncompressed_segments"], retained), name
        assert np.all(steps["compressed_segments"] >= 1), name
        assert np.all(steps["compressed_segments"] <= retained), name
        active = "compressed_segments" if run["compress"] else "uncompressed_segments"
        assert np.array_equal(steps["segments"], steps[active]), name
        if run.get("h_max", 0):
            assert run.get("h_max_reason"), (name, "A source-resolution cap needs a reason")
            assert max(steps["h"]) <= run["h_max"] + 1e-12, name
        error = np.r_[abs(samples["z"] - samples["exact_z"]),
                      abs(midpoints["z"] - midpoints["exact_z"])]
        if run.get("expected_missed", False):
            assert max(error) > 0.5, (name, float(max(error)))
        else:
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
            "peak_uncompressed_segments": int(max(steps["uncompressed_segments"])),
            "peak_compressed_segments": int(max(steps["compressed_segments"])),
        })
    with (HERE / "summary.csv").open("w") as stream:
        keys = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (HERE / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=REPO / "build")
    args = parser.parse_args()
    generate(args.build_dir.resolve())
    summarize()
    subprocess.run(["python3", "-B", str(HERE / "plot_gallery.py")], check=True)
