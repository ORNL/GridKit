"""Signal matching and validation metrics; no simulation or plotting."""

import numpy as np


def column_key(column):
    parts = column.split("_")
    body = parts[1:-1]
    if body and body[-1].lower() == parts[0].lower():
        body.pop()
    return " ".join(body)


def load_csv(path):
    with path.open() as stream:
        columns = stream.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[1] != len(columns) or not np.isfinite(data).all():
        raise ValueError(f"{path}: invalid or nonfinite signal data")
    return columns, data


def align(labels, values, expected):
    if (len(labels) != len(set(labels)) or len(expected) != len(set(expected))
            or sorted(labels) != sorted(expected)):
        raise ValueError("Output/reference channel mismatch")
    return values[:, [labels.index(label) for label in expected]]


def check_grid(time, values, reference_time, reference_values):
    if values.shape != reference_values.shape or time.shape != reference_time.shape:
        raise ValueError("Output/reference shape mismatch")
    if not np.allclose(time, reference_time, rtol=0.0, atol=1e-6):
        raise ValueError("Output/reference time-grid mismatch")
    if (not np.isfinite(values).all() or not np.isfinite(reference_values).all()
            or not np.isfinite(time).all() or not np.isfinite(reference_time).all()):
        raise ValueError("Nonfinite signals")
    if np.any(np.diff(reference_time) < 0) or reference_time[-1] <= reference_time[0]:
        raise ValueError("Invalid reference time grid")


def compare(columns, output, reference, suffix):
    expected, ref = load_csv(reference)
    indices = [i for i, column in enumerate(columns) if i > 0 and column.endswith(suffix)]
    values = align([column_key(columns[i]) for i in indices], output[:, indices], expected[1:])
    check_grid(output[:, 0], values, ref[:, 0], ref[:, 1:])
    error = values - ref[:, 1:]
    steps = np.diff(ref[:, 0])
    duration = ref[-1, 0] - ref[0, 0]
    # Temporal/channel RMSE, as used in the paper table (not the sweep's L2 norm).
    rmse = np.sqrt(np.median(steps[steps > 0]) / (values.shape[1] * duration)) * np.linalg.norm(error)
    scale = np.max(np.abs(ref[:, 1:]))
    if scale == 0:
        raise ValueError("Cannot normalize a zero reference signal")
    return (float(rmse), float(np.max(np.abs(error)) / scale)), expected, values
