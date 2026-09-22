"""Collect and read the optional segment-local IDA trace."""

import numpy as np

from common.simulation import disable_monitors, read_json, simulate, write_json


COUNTERS = ("accepted_steps", "residual_evals", "jacobian_evals", "error_test_failures")


def boundaries(study):
    return [0.0, *[event["time"] for event in study["events"]], study["tmax"]]


def collect(source, binary, work):
    study = read_json(source)
    case = read_json(source.parent / study["system_model_file"])
    disable_monitors(case)
    work = work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    write_json(work / "case.json", case)
    trace = (source.parent / study["solver_trace_file"]).resolve()
    pending = work / "trace.csv"
    pending.unlink(missing_ok=True)
    study.update(system_model_file="case.json", solver_trace_file=str(pending))
    _, elapsed = simulate(binary, work, study, "run")
    load(pending, boundaries(study))
    trace.parent.mkdir(parents=True, exist_ok=True)
    pending.replace(trace)
    return elapsed


def load(path, ends):
    trace = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                                       encoding="utf-8"))
    if trace.size == 0 or trace.dtype.names != ("segment", "phase", "t", "h", *COUNTERS):
        raise ValueError(f"{path}: invalid solver trace schema")
    if not np.isfinite(trace["t"]).all() or not np.isfinite(trace["h"]).all():
        raise ValueError(f"{path}: nonfinite times or step sizes")
    if set(trace["segment"]) != set(range(len(ends) - 1)):
        raise ValueError(f"{path}: incomplete event segments")
    for segment, (start, end) in enumerate(zip(ends, ends[1:])):
        rows = trace[trace["segment"] == segment]
        if (rows["phase"][0] != "init" or rows["phase"][-1] != "end"
                or rows["t"][0] != start or rows["t"][-1] != end
                or np.any(rows["phase"][1:-1] != "step")):
            raise ValueError(f"{path}: incomplete segment {segment}")
        steps = rows[1:-1]
        if (not steps.size or steps["t"][0] <= start or steps["t"][-1] < end
                or np.any(np.diff(steps["t"]) <= 0) or np.any(steps["h"] <= 0)
                or not np.array_equal(rows["accepted_steps"], np.r_[0, np.arange(1, len(steps)+1), len(steps)])):
            raise ValueError(f"{path}: invalid accepted steps in segment {segment}")
        for key in COUNTERS:
            if np.any(rows[key] < 0) or np.any(rows[key] != np.floor(rows[key])) or np.any(np.diff(rows[key]) < 0):
                raise ValueError(f"{path}: invalid {key} in segment {segment}")
    return trace


def step_segments(trace, ends):
    for segment, (start, end) in enumerate(zip(ends, ends[1:])):
        rows = trace[(trace["segment"] == segment) & (trace["phase"] == "step")]
        yield np.r_[start, np.minimum(rows["t"], end)], np.r_[rows["h"][0], rows["h"]]


def cumulative_counts(trace, ends):
    times, counts = [np.array([0.0])], [np.zeros((1, len(COUNTERS)), dtype=np.int64)]
    offset = np.zeros(len(COUNTERS), dtype=np.int64)
    for segment, end in enumerate(ends[1:]):
        rows = trace[trace["segment"] == segment]
        local = np.column_stack([rows[key] for key in COUNTERS])
        times.append(np.minimum(rows["t"], end))
        counts.append(local + offset)
        offset += local[-1]
    return np.concatenate(times), np.concatenate(counts)
