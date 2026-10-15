"""
SLURM workflow utilities for UQ experiments.

Provides reusable functions for:
1. SLURM script generation (sim, collect, resubmit jobs)
2. Job submission and tracking (sbatch calls, dependency handling)
3. Job status monitoring (sacct queries, parsing)
4. Failure analysis (parse failed run logs)
5. Metadata persistence (JSON I/O for kernel restart recovery)
6. Configuration utilities (slicing, formatting)

Designed to work with the notebook-based UQ workflow in uq_setup.ipynb.
Configuration (SLURM_ACCOUNT, WALLTIME, etc.) is passed as arguments,
keeping config visible and editable in the notebook.
"""

import os
import subprocess
import json
import glob
import datetime
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional, Union

# ============================================================================
# 6. Configuration & Utility Functions
# ============================================================================


def get_n_digits(n_samples: int) -> int:
    """
    Calculate number of digits needed for zero-padded run indices.

    Args:
        n_samples: Total number of samples

    Returns:
        Number of digits (e.g., 4 for n_samples=1025)
    """
    return len(str(n_samples - 1))


def compute_slices(
    n_samples: int,
    n_nodes: int,
) -> List[Tuple[int, int]]:
    """
    Divide samples into chunks for parallel SLURM nodes.

    Args:
        n_samples: Total number of samples
        n_nodes: Number of SLURM nodes to create

    Returns:
        List of (lo, hi) tuples, one per node. Each tuple is inclusive.
        Example: [(0, 255), (256, 511), ...] for n_samples=1025, n_nodes=4
    """
    chunk_size = (n_samples + n_nodes - 1) // n_nodes
    slices = [
        (i * chunk_size, min((i + 1) * chunk_size - 1, n_samples - 1))
        for i in range(n_nodes)
    ]
    return slices


def format_indices_str(indices: Union[List[int], range], n_digits: int) -> str:
    """
    Format list of indices as space-separated zero-padded string.

    Args:
        indices: List or range of integers
        n_digits: Width for zero-padding (e.g., 4 for "0042")

    Returns:
        Space-separated string (e.g., "0000 0001 0002 0003")
    """
    return " ".join(f"{i:0{n_digits}d}" for i in indices)


# ============================================================================
# 5. Metadata Persistence (JSON I/O)
# ============================================================================


def save_job_metadata(
    file_path: str,
    sim_job_ids: List[str],
    collect_job_id: Optional[str],
    uq_run_root: str,
    n_samples: int,
) -> None:
    """
    Save job metadata to JSON file for later recovery (kernel restart).

    Args:
        file_path: Path to save metadata (usually {UQ_RUN_ROOT}/slurm/job_ids.txt)
        sim_job_ids: List of simulation job IDs from sbatch submissions
        collect_job_id: Optional collection job ID (or None if not using COLLECT_IN_SLURM)
        uq_run_root: UQ experiment root directory
        n_samples: Total number of samples

    Raises:
        IOError: If file cannot be written
    """
    metadata = {
        "sim_job_ids": sim_job_ids,
        "collect_job_id": collect_job_id,
        "timestamp": datetime.datetime.now().isoformat(),
        "uq_run_root": uq_run_root,
        "n_samples": n_samples,
    }
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w") as f:
        json.dump(metadata, f, indent=2)


def load_job_metadata(file_path: str) -> Dict:
    """
    Load job metadata from JSON file.

    Args:
        file_path: Path to metadata file

    Returns:
        Dictionary with keys: sim_job_ids, collect_job_id, timestamp, uq_run_root, n_samples

    Raises:
        FileNotFoundError: If file does not exist
        json.JSONDecodeError: If file is not valid JSON
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Job metadata file not found: {file_path}")

    with open(file_path, "r") as f:
        metadata = json.load(f)

    return metadata


def save_resubmit_selection(
    file_path: str,
    indices: List[int],
) -> None:
    """
    Save list of run indices marked for resubmission.

    Args:
        file_path: Path to save indices (usually {UQ_RUN_ROOT}/slurm/resubmit_indices.txt)
        indices: List of run indices to resubmit
    """
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w") as f:
        f.write("\n".join(str(i) for i in indices))


def load_resubmit_selection(file_path: str) -> List[int]:
    """
    Load list of run indices marked for resubmission.

    Args:
        file_path: Path to resubmit indices file

    Returns:
        List of integers

    Raises:
        FileNotFoundError: If file does not exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Resubmit indices file not found: {file_path}")

    with open(file_path, "r") as f:
        indices = [int(line.strip()) for line in f if line.strip()]

    return indices


def save_resubmit_metadata(
    file_path: str,
    resubmit_job_ids: List[str],
    resubmit_indices: List[int],
) -> None:
    """
    Save resubmit job metadata to JSON file.

    Args:
        file_path: Path to save metadata (usually {UQ_RUN_ROOT}/slurm/resubmit_job_ids.txt)
        resubmit_job_ids: List of resubmit job IDs
        resubmit_indices: List of run indices being resubmitted
    """
    metadata = {
        "resubmit_job_ids": resubmit_job_ids,
        "resubmit_indices": resubmit_indices,
        "timestamp": datetime.datetime.now().isoformat(),
    }
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w") as f:
        json.dump(metadata, f, indent=2)


def load_resubmit_metadata(file_path: str) -> Dict:
    """
    Load resubmit job metadata from JSON file.

    Args:
        file_path: Path to resubmit metadata file

    Returns:
        Dictionary with keys: resubmit_job_ids, resubmit_indices, timestamp

    Raises:
        FileNotFoundError: If file does not exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Resubmit metadata file not found: {file_path}")

    with open(file_path, "r") as f:
        metadata = json.load(f)

    return metadata


# ============================================================================
# 4. Failure Analysis
# ============================================================================


def parse_failed_runs(
    uq_run_root: str,
    n_samples: int,
) -> Tuple[Set[int], Dict[str, List[Tuple[int, str]]]]:
    """
    Scan chunk_*_failed.txt files and extract failed run indices.

    Args:
        uq_run_root: UQ experiment root directory
        n_samples: Total number of samples (used for context only)

    Returns:
        Tuple of:
            - Set of failed run indices (e.g., {5, 17, 42})
            - Dict mapping chunk_name -> list of (run_num, status) tuples
              (e.g., {"00": [(5, "FAILED"), (17, "NO_SOLVER")], ...})
    """
    slurm_dir = os.path.join(uq_run_root, "slurm")
    failed_runs = set()
    failed_by_chunk = {}

    for chunk_file in sorted(glob.glob(os.path.join(slurm_dir, "chunk_*_failed.txt"))):
        chunk_name = (
            os.path.basename(chunk_file)
            .replace("chunk_", "")
            .replace("_failed.txt", "")
        )
        failed_by_chunk[chunk_name] = []

        try:
            with open(chunk_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if line:
                        parts = line.split()
                        if len(parts) >= 2:
                            status = parts[0]
                            run_dir = parts[1]
                            # Extract run number from "sample_runs/run_NNNN"
                            if "run_" in run_dir:
                                run_num_str = run_dir.split("run_")[-1].rstrip("/")
                                try:
                                    run_num = int(run_num_str)
                                    failed_runs.add(run_num)
                                    failed_by_chunk[chunk_name].append(
                                        (run_num, status)
                                    )
                                except ValueError:
                                    pass
        except FileNotFoundError:
            pass

    return failed_runs, failed_by_chunk


def print_failure_summary(
    failed_runs: Set[int],
    failed_by_chunk: Dict[str, List[Tuple[int, str]]],
    n_samples: int,
) -> None:
    """
    Pretty-print failure summary to console.

    Args:
        failed_runs: Set of failed run indices
        failed_by_chunk: Dict mapping chunk_name -> [(run_num, status), ...]
        n_samples: Total number of samples (for formatting)
    """
    n_digits = get_n_digits(n_samples)

    if not failed_runs:
        print(f"✓ No failures found! All runs completed successfully.")
    else:
        print(
            f"✗ Found {len(failed_runs)} failed runs across {len(failed_by_chunk)} chunks:"
        )
        print()
        for chunk_name in sorted(failed_by_chunk.keys()):
            runs = failed_by_chunk[chunk_name]
            if runs:
                print(f"  chunk_{chunk_name}: {len(runs)} failures")
                for run_num, status in sorted(runs)[:5]:  # Show first 5
                    print(f"    - run_{run_num:0{n_digits}d}: {status}")
                if len(runs) > 5:
                    print(f"    ... and {len(runs) - 5} more")
        print()
        print(f"Failed run indices: {sorted(failed_runs)}")


# ============================================================================
# 3. Job Status Monitoring (sacct queries)
# ============================================================================


def parse_sacct_output(stdout_text: str) -> Dict[str, int]:
    """
    Parse sacct output and count job states (by unique base job ID).

    Note: sacct returns 3 entries per job (main, .batch, .extern). This function
    extracts the base job ID and counts states only once per unique job.

    Args:
        stdout_text: Raw stdout from `sacct -j ...` command

    Returns:
        Dictionary with keys: "completed", "running", "pending", "failed"
        Example: {"completed": 4, "running": 0, "pending": 0, "failed": 0}
    """
    counts = {"completed": 0, "running": 0, "pending": 0, "failed": 0}

    # Track unique base job IDs (strip .batch, .extern) and their states
    # Use a dict to store the best state for each base job ID
    job_states = {}

    lines = stdout_text.strip().split("\n")[1:]  # Skip header
    for line in lines:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 3:
            job_id_str = parts[0]
            state = parts[2]

            # Extract base job ID (e.g., "12345" from "12345.batch" or "12345.extern")
            base_id = job_id_str.split(".")[0]

            # Keep the "best" state (prioritize COMPLETED over others for the base ID)
            # This way, if we see both 12345.batch and 12345.extern, we count once
            if base_id not in job_states:
                job_states[base_id] = state
            else:
                # Prefer COMPLETED if we see it on any variant
                if state == "COMPLETED" or (
                    job_states[base_id] != "COMPLETED" and state != "COMPLETED"
                ):
                    job_states[base_id] = state

    # Now count unique states
    for base_id, state in job_states.items():
        if state == "COMPLETED":
            counts["completed"] += 1
        elif state in ("FAILED", "TIMEOUT", "CANCELLED"):
            counts["failed"] += 1
        elif state == "RUNNING":
            counts["running"] += 1
        elif state in ("PENDING", "CONFIGURING"):
            counts["pending"] += 1

    return counts


def check_job_status(
    job_ids: List[str],
    collect_job_id: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, Union[Dict[str, int], str]]:
    """
    Query sacct for job status and print summary.

    Args:
        job_ids: List of job IDs (e.g., ["12345", "12346"])
        collect_job_id: Optional collection job ID (will query separately)
        verbose: If True, print sacct output; if False, only return counts

    Returns:
        Dictionary with keys:
            - "sim_counts": {"completed", "running", "pending", "failed"}
            - "collect_state": "COMPLETED" or None if no collect job
    """
    result = {"sim_counts": None, "collect_state": None}

    if not job_ids:
        print("No job IDs to check")
        return result

    # Query sim jobs
    cmd = [
        "sacct",
        "-j",
        ",".join(job_ids),
        "--format=JobID,JobName,State,Elapsed,ExitCode",
    ]
    if verbose:
        print(" ".join(cmd))

    proc_result = subprocess.run(cmd, capture_output=True, text=True)
    if proc_result.returncode == 0:
        if verbose:
            print(proc_result.stdout)
        counts = parse_sacct_output(proc_result.stdout)
        result["sim_counts"] = counts

        if verbose:
            total_jobs = len(job_ids)
            print(
                f"\nSim Jobs Status ({counts['completed']}/{total_jobs} COMPLETED): "
                f"{counts['running']} RUNNING, {counts['pending']} PENDING, "
                f"{counts['failed']} FAILED"
            )
    else:
        print(f"ERROR querying jobs: {proc_result.stderr}")

    # Query collect job if provided
    if collect_job_id:
        if verbose:
            print(f"\n--- Collect Job ---")
        cmd = [
            "sacct",
            "-j",
            collect_job_id,
            "--format=JobID,JobName,State,Elapsed,ExitCode",
        ]
        proc_result = subprocess.run(cmd, capture_output=True, text=True)
        if verbose:
            print(proc_result.stdout)
        # Extract state from first line of output (after header)
        lines = proc_result.stdout.strip().split("\n")
        if len(lines) > 1:
            parts = lines[1].split()
            if len(parts) >= 3:
                state = parts[2]
                result["collect_state"] = state
                if verbose:
                    print(f"Collect Job Status: {state}")

    return result


def get_experiment_status(
    sim_job_ids: List[str],
    collect_job_id: Optional[str] = None,
) -> str:
    """
    Query sacct to determine overall experiment status at a glance.

    Useful for multi-experiment dashboards where you need a single status string
    per experiment (e.g., "✅ COMPLETED", "🔄 IN PROGRESS", "❌ FAILED", etc.).

    Args:
        sim_job_ids: List of simulation job IDs
        collect_job_id: Optional collection job ID

    Returns:
        Status string: "✅ COMPLETED", "🔄 IN PROGRESS", "❌ FAILED",
                       "⚠ CANCELLED", "⏳ PENDING", or "⚠ UNKNOWN"
    """
    all_job_ids = sim_job_ids.copy()
    if collect_job_id is not None:
        all_job_ids.append(collect_job_id)

    # Query sacct for all job states
    try:
        job_ids_str = ",".join(str(jid) for jid in all_job_ids)
        result = subprocess.run(
            [
                "sacct",
                "-j",
                job_ids_str,
                "--format=JobID,State",
                "--parsable2",
                "--noheader",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode != 0:
            return "⚠ UNKNOWN"

        # Parse states from sacct output
        states = {}
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("|")
            if len(parts) >= 2:
                job_id_str, state = parts[0].strip(), parts[1].strip()
                # Extract base job ID (strip .batch, .extern suffixes)
                base_id = job_id_str.split(".")[0]
                states[base_id] = state

        if not states:
            return "⏳ SUBMITTED"

        # Determine overall status
        has_running = any("RUNNING" in s for s in states.values())
        has_failed = any("FAILED" in s for s in states.values())
        has_cancelled = any("CANCELLED" in s for s in states.values())
        all_completed = all("COMPLETED" in s for s in states.values())

        if has_failed:
            return "❌ FAILED"
        elif has_cancelled:
            return "⚠ CANCELLED"
        elif has_running:
            return "🔄 IN PROGRESS"
        elif all_completed:
            return "✅ COMPLETED"
        else:
            return "⏳ PENDING"

    except Exception:
        return "⚠ UNKNOWN"


# ============================================================================
# 2. Job Submission & Tracking
# ============================================================================


def submit_job(script_path: str) -> Optional[str]:
    """
    Submit a single sbatch script and return job ID.

    Args:
        script_path: Path to .sh script file

    Returns:
        Job ID (e.g., "12345") or None if submission failed
    """
    cmd = ["sbatch", script_path]
    print(" ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
        return None

    # Extract job ID from stdout (e.g., "Submitted batch job 12345")
    jid = result.stdout.strip().split()[-1]
    print(f"  submitted job {jid}")
    return jid


def submit_jobs_batch(script_paths: List[str]) -> List[str]:
    """
    Submit multiple sbatch scripts and collect job IDs.

    Args:
        script_paths: List of paths to .sh script files

    Returns:
        List of job IDs (empty if all submissions failed)
    """
    job_ids = []
    for script_path in script_paths:
        jid = submit_job(script_path)
        if jid:
            job_ids.append(jid)

    print(f"\nSubmitted {len(job_ids)} jobs: {job_ids}")
    return job_ids


def submit_collect_with_dependency(
    collect_script_path: str,
    sim_job_ids: List[str],
) -> Optional[str]:
    """
    Submit collect job with dependency on sim jobs (--dependency=afterok:...).

    Args:
        collect_script_path: Path to collect.sh script
        sim_job_ids: List of sim job IDs to depend on

    Returns:
        Collect job ID or None if submission failed

    Raises:
        ValueError: If sim_job_ids is empty
    """
    if not sim_job_ids:
        raise ValueError("Cannot submit collect job without sim job IDs")

    dep = "afterok:" + ":".join(sim_job_ids)
    cmd = ["sbatch", f"--dependency={dep}", collect_script_path]
    print(f"\n{' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR submitting collect job: {result.stderr.strip()}")
        return None

    collect_jid = result.stdout.strip().split()[-1]
    print(f"  submitted collect job {collect_jid}  (runs after all sim jobs complete)")
    return collect_jid


# ============================================================================
# 1. SLURM Script Generation (Template Abstraction)
# ============================================================================


def generate_sim_script(
    node_i: int,
    indices: List[int],
    uq_run_root: str,
    slurm_dir: str,
    gridkit_build_dir: str,
    gridkit_py_utils: str,
    gridkit_runner: str,
    workers_per_node: int,
    walltime: str,
    slurm_account: str,
    slurm_partition: str,
    slurm_qos: Optional[str] = None,
    slurm_modules: Optional[List[str]] = None,
    n_digits: int = 4,
) -> str:
    """
    Generate bash script for simulation job on one node.

    Args:
        node_i: Node index (0, 1, 2, ...)
        indices: List of run indices to process on this node (already formatted as ints)
        uq_run_root: UQ experiment root directory
        slurm_dir: Directory to save logs/failed.txt files
        gridkit_build_dir: GRIDKIT build directory (contains bin/)
        gridkit_py_utils: GRIDKIT py-utils directory (for gridkit_utils import)
        gridkit_runner: DynamicSimulation command (full path or just name if in PATH)
        workers_per_node: Number of parallel workers (-P flag to xargs)
        walltime: Wall time limit (e.g., "01:00:00")
        slurm_account: Account name (--account)
        slurm_partition: Partition name (--partition)
        slurm_qos: Optional QOS (--qos)
        slurm_modules: Optional list of modules to load (e.g., ["intel-mkl", "gcc"])
        n_digits: Number of digits for zero-padding run indices

    Returns:
        Bash script as string
    """
    indices_str = format_indices_str(indices, n_digits)

    module_loads = ""
    if slurm_modules:
        module_loads = "module load " + " ".join(slurm_modules)

    qos_line = ""
    if slurm_qos:
        qos_line = f"#SBATCH --qos={slurm_qos}"

    script = f"""#!/bin/bash
#SBATCH --job-name="sim_{node_i:02d}"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={workers_per_node}
#SBATCH --time={walltime}
#SBATCH --account={slurm_account}
#SBATCH --partition={slurm_partition}
#SBATCH --output={slurm_dir}/chunk_{node_i:02d}.log
{qos_line}

set -e
export PATH={gridkit_build_dir}/bin:$PATH
export PYTHONPATH={gridkit_py_utils}:$PYTHONPATH

# Module loads (if any)
{module_loads}

cd {uq_run_root}

INDICES=({indices_str})
printf '%s\\n' "${{INDICES[@]}}" | xargs -P {workers_per_node} -I {{}} bash -c '
    dir="sample_runs/run_{{}}"
    solver=$(ls "$dir"/*.solver.json 2>/dev/null | head -1)
    if [[ -z "$solver" ]]; then
        echo "NO_SOLVER $dir" >> "{slurm_dir}/chunk_{node_i:02d}_failed.txt"
        exit 0
    fi
    cd "$dir"
    {gridkit_runner} "$(basename $solver)" > stdout.txt 2> stderr.txt \\
        || echo "FAILED $dir" >> "{slurm_dir}/chunk_{node_i:02d}_failed.txt"
'
"""
    return script


def generate_collect_script(
    uq_run_root: str,
    uq_out_path: str,
    slurm_dir: str,
    gridkit_py_utils: str,
    workers: int,
    walltime: str,
    slurm_account: str,
    slurm_partition: str,
    conda_env: str,
    n_samples: int,
    slurm_qos: Optional[str] = None,
) -> str:
    """
    Generate bash script for collection job (aggregates results after sims complete).

    Args:
        uq_run_root: UQ experiment root directory
        uq_out_path: Output path for aggregated results
        slurm_dir: Directory for logs
        gridkit_py_utils: GRIDKIT py-utils directory
        workers: Number of parallel workers for collection
        walltime: Wall time limit
        slurm_account: Account name
        slurm_partition: Partition name
        conda_env: Path to conda environment
        n_samples: Total number of samples
        slurm_qos: Optional QOS

    Returns:
        Bash script as string
    """
    qos_line = ""
    if slurm_qos:
        qos_line = f"#SBATCH --qos={slurm_qos}"

    script = f"""#!/bin/bash
#SBATCH --job-name="collect"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={workers}
#SBATCH --time={walltime}
#SBATCH --account={slurm_account}
#SBATCH --partition={slurm_partition}
#SBATCH --output={slurm_dir}/collect.log
{qos_line}

source /nopt/nrel/apps/software/mamba/etc/profile.d/conda.sh
conda activate {conda_env}

cd {uq_run_root}
python3 << 'PYEOF'
import sys
sys.path.insert(0, "{gridkit_py_utils}")
from gridkit_utils import collect_parallel

result, missing = collect_parallel(
    "{uq_run_root}",
    "{uq_out_path}",
    n_workers={workers},
    n_samples={n_samples},
)
print(f"Collected {{len(result)}} results")
if missing:
    print(f"Missing runs: {{missing}}")
    sys.exit(1)
PYEOF
"""
    return script


def generate_resubmit_script(
    node_i: int,
    indices: List[int],
    uq_run_root: str,
    slurm_dir: str,
    gridkit_build_dir: str,
    gridkit_py_utils: str,
    gridkit_runner: str,
    workers_per_node: int,
    walltime: str,
    slurm_account: str,
    slurm_partition: str,
    slurm_qos: Optional[str] = None,
    slurm_modules: Optional[List[str]] = None,
    n_digits: int = 4,
) -> str:
    """
    Generate bash script for resubmission job (failed runs only).

    Similar to generate_sim_script but for a subset of indices.

    Args:
        node_i: Node index for resubmit chunk
        indices: List of run indices to resubmit
        uq_run_root: UQ experiment root directory
        slurm_dir: Directory for logs
        gridkit_build_dir: GRIDKIT build directory
        gridkit_py_utils: GRIDKIT py-utils directory
        gridkit_runner: DynamicSimulation command
        workers_per_node: Number of parallel workers
        walltime: Wall time limit
        slurm_account: Account name
        slurm_partition: Partition name
        slurm_qos: Optional QOS
        slurm_modules: Optional list of modules to load
        n_digits: Number of digits for zero-padding

    Returns:
        Bash script as string
    """
    indices_str = format_indices_str(indices, n_digits)

    module_loads = ""
    if slurm_modules:
        module_loads = "module load " + " ".join(slurm_modules)

    qos_line = ""
    if slurm_qos:
        qos_line = f"#SBATCH --qos={slurm_qos}"

    script = f"""#!/bin/bash
#SBATCH --job-name="resub_{node_i:02d}"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={workers_per_node}
#SBATCH --time={walltime}
#SBATCH --account={slurm_account}
#SBATCH --partition={slurm_partition}
#SBATCH --output={slurm_dir}/resubmit_{node_i:02d}.log
{qos_line}

set -e
export PATH={gridkit_build_dir}/bin:$PATH
export PYTHONPATH={gridkit_py_utils}:$PYTHONPATH

# Module loads (if any)
{module_loads}

cd {uq_run_root}

INDICES=({indices_str})
printf '%s\\n' "${{INDICES[@]}}" | xargs -P {workers_per_node} -I {{}} bash -c '
    dir="sample_runs/run_{{}}"
    solver=$(ls "$dir"/*.solver.json 2>/dev/null | head -1)
    if [[ -z "$solver" ]]; then
        echo "NO_SOLVER $dir" >> "{slurm_dir}/resubmit_{node_i:02d}_failed.txt"
        exit 0
    fi
    cd "$dir"
    {gridkit_runner} "$(basename $solver)" > stdout.txt 2> stderr.txt \\
        || echo "FAILED $dir" >> "{slurm_dir}/resubmit_{node_i:02d}_failed.txt"
'
"""
    return script


def write_script_files(
    uq_run_root: str,
    slices: List[Tuple[int, int]],
    slurm_account: str,
    slurm_partition: str,
    walltime: str,
    workers_per_node: int,
    gridkit_build_dir: str,
    gridkit_py_utils: str,
    gridkit_runner: str,
    gridkit_runner_name: str,
    n_samples: int,
    slurm_qos: Optional[str] = None,
    slurm_modules: Optional[List[str]] = None,
    collect_in_slurm: bool = False,
    collect_workers: int = 64,
    collect_walltime: str = "02:00:00",
    conda_env: Optional[str] = None,
    uq_out_path: Optional[str] = None,
) -> Tuple[List[str], Optional[str]]:
    """
    Generate and write all sbatch scripts to disk.

    Args:
        uq_run_root: UQ experiment root directory
        slices: List of (lo, hi) tuples from compute_slices()
        slurm_account: --account value
        slurm_partition: --partition value
        walltime: --time value (e.g., "01:00:00")
        workers_per_node: -P value for xargs (threads per node)
        gridkit_build_dir: Path to GRIDKIT build
        gridkit_py_utils: Path to py-utils (for PYTHONPATH)
        gridkit_runner: Full path to DynamicSimulation (or just name if in PATH)
        gridkit_runner_name: Display name of runner for docs
        n_samples: Total number of samples (for padding calculation)
        slurm_qos: Optional --qos value
        slurm_modules: Optional list of modules to load
        collect_in_slurm: If True, also generate collect.sh
        collect_workers: Threads for collection job
        collect_walltime: Wall time for collection job
        conda_env: Path to conda environment (for collect job)
        uq_out_path: Output path for results (for collect job)

    Returns:
        Tuple of:
            - List of sim script paths
            - Collect script path (or None if not collect_in_slurm)
    """
    n_digits = get_n_digits(n_samples)
    slurm_dir = os.path.join(uq_run_root, "slurm")
    os.makedirs(slurm_dir, exist_ok=True)

    script_paths = []

    # Generate sim scripts
    for node_i, (lo, hi) in enumerate(slices):
        indices = list(range(lo, hi + 1))
        script_content = generate_sim_script(
            node_i=node_i,
            indices=indices,
            uq_run_root=uq_run_root,
            slurm_dir=slurm_dir,
            gridkit_build_dir=gridkit_build_dir,
            gridkit_py_utils=gridkit_py_utils,
            gridkit_runner=gridkit_runner,
            workers_per_node=workers_per_node,
            walltime=walltime,
            slurm_account=slurm_account,
            slurm_partition=slurm_partition,
            slurm_qos=slurm_qos,
            slurm_modules=slurm_modules,
            n_digits=n_digits,
        )

        script_path = os.path.join(slurm_dir, f"chunk_{node_i:02d}.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        script_paths.append(script_path)
        print(f"  wrote {script_path}  (runs {lo:04d}-{hi:04d})")

    # Generate collect script if requested
    collect_script_path = None
    if collect_in_slurm:
        collect_content = generate_collect_script(
            uq_run_root=uq_run_root,
            uq_out_path=uq_out_path or os.path.join(uq_run_root, "results"),
            slurm_dir=slurm_dir,
            gridkit_py_utils=gridkit_py_utils,
            workers=collect_workers,
            walltime=collect_walltime,
            slurm_account=slurm_account,
            slurm_partition=slurm_partition,
            conda_env=conda_env or "",
            n_samples=n_samples,
            slurm_qos=slurm_qos,
        )

        collect_script_path = os.path.join(slurm_dir, "collect.sh")
        with open(collect_script_path, "w") as f:
            f.write(collect_content)
        os.chmod(collect_script_path, 0o755)
        print(
            f"    wrote {collect_script_path}  (collect job, {collect_workers} threads)"
        )

    print(
        f"\nTotal nodes: {len(script_paths)} sim"
        + (f" + 1 collect = {len(script_paths) + 1}" if collect_in_slurm else "")
    )
    print(f"\nReview sim script:     cat {slurm_dir}/chunk_00.sh")
    if collect_script_path:
        print(f"Review collect script: cat {collect_script_path}")

    return script_paths, collect_script_path


def generate_and_write_resubmit_scripts(
    uq_run_root: str,
    resubmit_indices: List[int],
    n_nodes: int,
    slurm_account: str,
    slurm_partition: str,
    walltime: str,
    workers_per_node: int,
    gridkit_build_dir: str,
    gridkit_py_utils: str,
    gridkit_runner: str,
    n_samples: int,
    slurm_qos: Optional[str] = None,
    slurm_modules: Optional[List[str]] = None,
) -> List[str]:
    """
    Generate and write resubmit scripts for failed runs.

    Args:
        uq_run_root: UQ experiment root directory
        resubmit_indices: List of run indices to resubmit
        n_nodes: Number of resubmit nodes to create
        slurm_account: --account value
        slurm_partition: --partition value
        walltime: --time value
        workers_per_node: -P value for xargs
        gridkit_build_dir: Path to GRIDKIT build
        gridkit_py_utils: Path to py-utils
        gridkit_runner: Full path to DynamicSimulation
        n_samples: Total number of samples (for padding)
        slurm_qos: Optional --qos value
        slurm_modules: Optional list of modules to load

    Returns:
        List of resubmit script paths
    """
    n_digits = get_n_digits(n_samples)
    slurm_dir = os.path.join(uq_run_root, "slurm")
    os.makedirs(slurm_dir, exist_ok=True)

    # Divide resubmit indices into chunks (like original slices)
    chunk_size = (len(resubmit_indices) + n_nodes - 1) // n_nodes

    script_paths = []

    for node_i in range(n_nodes):
        lo_idx = node_i * chunk_size
        hi_idx = min((node_i + 1) * chunk_size, len(resubmit_indices))

        if lo_idx >= hi_idx:
            break

        chunk_indices = resubmit_indices[lo_idx:hi_idx]

        script_content = generate_resubmit_script(
            node_i=node_i,
            indices=chunk_indices,
            uq_run_root=uq_run_root,
            slurm_dir=slurm_dir,
            gridkit_build_dir=gridkit_build_dir,
            gridkit_py_utils=gridkit_py_utils,
            gridkit_runner=gridkit_runner,
            workers_per_node=workers_per_node,
            walltime=walltime,
            slurm_account=slurm_account,
            slurm_partition=slurm_partition,
            slurm_qos=slurm_qos,
            slurm_modules=slurm_modules,
            n_digits=n_digits,
        )

        script_path = os.path.join(slurm_dir, f"resubmit_{node_i:02d}.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        script_paths.append(script_path)
        print(f"  wrote {script_path}  ({len(chunk_indices)} runs)")

    return script_paths
