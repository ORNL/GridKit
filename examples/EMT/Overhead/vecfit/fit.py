from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXAMPLE_DIR = HERE.parent
INPUT_CSV = EXAMPLE_DIR / "overhead.response.csv"
OUTPUT_DIR = EXAMPLE_DIR / "output"
YC_CSV = OUTPUT_DIR / "yc.csv"
HMIN_CSV = OUTPUT_DIR / "hmin.csv"
DELAY_JSON = OUTPUT_DIR / "delay.json"
YC_MODEL_JSON = OUTPUT_DIR / "yc.model.json"
HMIN_MODEL_JSON = OUTPUT_DIR / "hmin.model.json"
YC_POLES = 10
HMIN_POLES = 20

YC_COLUMN = re.compile(r"^Overhead_Yc_(real|imag)_(\d+)_(\d+)$")


def coordinate_column(fieldnames: list[str]) -> str:
    if "omega" in fieldnames:
        return "omega"
    if "t" in fieldnames:
        return "t"
    raise ValueError(f"{INPUT_CSV} must contain an omega column")


def read_monitor_csv() -> tuple[list[str], str, list[dict[str, str]]]:
    with INPUT_CSV.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{INPUT_CSV} must contain a header")

        fieldnames = reader.fieldnames
        omega_column = coordinate_column(fieldnames)
        return fieldnames, omega_column, list(reader)


def yc_columns(fieldnames: list[str]) -> tuple[int, int, list[tuple[int, int]]]:
    columns: dict[str, set[tuple[int, int]]] = {"real": set(), "imag": set()}

    for name in fieldnames:
        match = YC_COLUMN.match(name)
        if match is None:
            continue
        part, row, col = match.groups()
        columns[part].add((int(row), int(col)))

    real = columns["real"]
    imag = columns["imag"]
    if not real:
        raise ValueError(f"{INPUT_CSV} contains no Overhead_Yc monitor columns")
    if real != imag:
        missing_imag = sorted(real - imag)
        missing_real = sorted(imag - real)
        raise ValueError(
            f"Yc real/imag columns are mismatched: "
            f"missing imag {missing_imag}, missing real {missing_real}"
        )

    rows = max(row for row, _ in real) + 1
    cols = max(col for _, col in real) + 1
    expected = {(row, col) for row in range(rows) for col in range(cols)}
    if real != expected:
        missing = sorted(expected - real)
        raise ValueError(f"Yc matrix columns are not complete: missing {missing}")

    return rows, cols, sorted(real)


def h_modes(fieldnames: list[str]) -> list[int]:
    real_prefix = "Overhead_H_real_"
    imag_prefix = "Overhead_H_imag_"
    tau_prefix = "Overhead_Tau_"

    real = {
        int(name.removeprefix(real_prefix))
        for name in fieldnames
        if name.startswith(real_prefix)
    }
    imag = {
        int(name.removeprefix(imag_prefix))
        for name in fieldnames
        if name.startswith(imag_prefix)
    }
    tau = {
        int(name.removeprefix(tau_prefix))
        for name in fieldnames
        if name.startswith(tau_prefix)
    }

    if not real:
        raise ValueError(f"{INPUT_CSV} contains no Overhead_H monitor columns")
    if real != imag or real != tau:
        raise ValueError(
            f"H real/imag/Tau columns are mismatched: "
            f"missing imag {sorted(real - imag)}, missing real {sorted(imag - real)}, "
            f"missing tau {sorted(real - tau)}"
        )

    return sorted(real)


def frequency_hz(input_row: dict[str, str], omega_column: str) -> float:
    return float(input_row[omega_column]) / (2.0 * math.pi)


def hmin_value(input_row: dict[str, str], mode: int, tau_min: float, omega: float) -> complex:
    real = float(input_row[f"Overhead_H_real_{mode}"])
    imag = float(input_row[f"Overhead_H_imag_{mode}"])
    angle = omega * tau_min
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    return complex(
        real * cos_angle - imag * sin_angle,
        real * sin_angle + imag * cos_angle,
    )


def write_vecfit_yc_csv(
    monitor_rows: list[dict[str, str]],
    omega_column: str,
    entries: list[tuple[int, int]],
) -> None:
    with YC_CSV.open("w", newline="") as out_stream:
        writer = csv.writer(out_stream)
        header = ["freq_Hz"]
        for row, col in entries:
            header.extend([f"re_Yc_{row}_{col}", f"im_Yc_{row}_{col}"])
        writer.writerow(header)

        for input_row in monitor_rows:
            output_row = [f"{frequency_hz(input_row, omega_column):.17e}"]
            for row, col in entries:
                output_row.append(input_row[f"Overhead_Yc_real_{row}_{col}"])
                output_row.append(input_row[f"Overhead_Yc_imag_{row}_{col}"])
            writer.writerow(output_row)


def write_vecfit_hmin_csv(
    monitor_rows: list[dict[str, str]],
    omega_column: str,
    modes: list[int],
) -> list[float]:
    tau_min_by_mode = {
        mode: min(float(row[f"Overhead_Tau_{mode}"]) for row in monitor_rows)
        for mode in modes
    }

    with HMIN_CSV.open("w", newline="") as out_stream:
        writer = csv.writer(out_stream)
        header = ["freq_Hz"]
        for mode in modes:
            header.extend([f"re_Hmin_{mode}_0", f"im_Hmin_{mode}_0"])
        writer.writerow(header)

        for input_row in monitor_rows:
            omega = float(input_row[omega_column])
            output_row = [f"{omega / (2.0 * math.pi):.17e}"]
            for mode in modes:
                value = hmin_value(input_row, mode, tau_min_by_mode[mode], omega)
                output_row.extend([f"{value.real:.17e}", f"{value.imag:.17e}"])
            writer.writerow(output_row)

    return [tau_min_by_mode[mode] for mode in modes]


def write_vecfit_inputs() -> None:
    fieldnames, omega_column, monitor_rows = read_monitor_csv()
    yc_rows, yc_cols, yc_entries = yc_columns(fieldnames)
    modes = h_modes(fieldnames)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_vecfit_yc_csv(monitor_rows, omega_column, yc_entries)
    tau_min = write_vecfit_hmin_csv(monitor_rows, omega_column, modes)
    DELAY_JSON.write_text(json.dumps(tau_min, indent=2) + "\n")

    print(f"wrote {yc_rows}x{yc_cols} Yc vecfit CSV: {YC_CSV}")
    print(f"wrote {len(modes)}x1 Hmin vecfit CSV: {HMIN_CSV}")
    print(f"wrote modal tau_min delays: {DELAY_JSON}")

    fit_vecfit_model(YC_CSV, YC_MODEL_JSON, yc_rows, yc_cols, YC_POLES)
    fit_vecfit_model(
        HMIN_CSV,
        HMIN_MODEL_JSON,
        len(modes),
        1,
        HMIN_POLES,
        terms="none",
        state_space=True,
    )


def fit_vecfit_model(
    csv_path: Path,
    model_path: Path,
    rows: int,
    cols: int,
    poles: int,
    terms: str | None = None,
    state_space: bool = False,
) -> None:
    vecfit = shutil.which("vecfit")
    if vecfit is None:
        cargo_vecfit = Path.home() / ".cargo" / "bin" / "vecfit"
        if cargo_vecfit.exists():
            vecfit = str(cargo_vecfit)
    if vecfit is None:
        raise RuntimeError("vecfit CLI not found on PATH")

    command = [
        vecfit,
        "fit",
        str(csv_path),
        "--matrix",
        str(rows),
        str(cols),
        "--poles",
        str(poles),
        "--weighted",
    ]
    if terms is not None:
        command.extend(["--terms", terms])
    if state_space:
        command.append("--state-space")
    command.extend(["--output", str(model_path)])

    subprocess.run(command, check=True)
    print(f"wrote {rows}x{cols} {poles}-pole vecfit model: {model_path}")


write_vecfit_inputs()
