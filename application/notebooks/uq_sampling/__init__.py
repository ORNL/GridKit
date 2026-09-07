"""Parameter sampling, run orchestration, and dispatch-patching for UQ sweeps.

Not yet wired into a notebook in this repository — added ahead of the UQ
sweep notebook that will exercise it. See application/notebooks/README.md.
"""

from .dispatch import patch_genrou_dispatch, plot_genrou_dispatch, read_genrou_dispatch
from .runs import MONITORABLE_VARS_BY_ELEMENT, collect_and_save, make_run_dir, run_sample
from .sampling import generate_samples

__all__ = [
    "MONITORABLE_VARS_BY_ELEMENT",
    "collect_and_save",
    "generate_samples",
    "make_run_dir",
    "patch_genrou_dispatch",
    "plot_genrou_dispatch",
    "read_genrou_dispatch",
    "run_sample",
]
