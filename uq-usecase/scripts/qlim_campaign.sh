#!/bin/bash
# Run a full campaign design (baseline + sited rows) with PV->PQ Q-limit
# switching enforced in the AC PF, at a given dc alphaP. Unlike highp_rerun.sh,
# baseline rows are kept: with Q limits the no-datacenter case also changes.
# usage: sbatch qlim_campaign.sh <src-campaign> <dst-campaign> <alphaP> <workers> [switch|acopf|acopf_noq]
#   e.g. sbatch qlim_campaign.sh faultscreen1 faultscreen1_ap08_qlim 0.8 16
#        sbatch qlim_campaign.sh faultscreen1 faultscreen1_ap08_acopf 0.8 16 acopf
# switch (default): DC OPF + AC PF with PV->PQ switching. acopf: PM.jl AC OPF.
# acopf_noq: the same AC OPF with generator Q limits widened (control).
#SBATCH --job-name=qlim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short

set -u
SRC=$1; DST=$2; ALPHA=$3; W=$4; METHOD=${5:-switch}
case "$METHOD" in
  switch) FLAGS=(--enforce-q-lims) ;;
  acopf)  FLAGS=(--ac-opf) ;;
  acopf_noq) FLAGS=(--ac-opf --ac-opf-relax-q) ;;
  *) echo "unknown method $METHOD"; exit 2 ;;
esac
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
ROOT=/kfs2/projects/scidac/scidac-data/gridkit-campaigns
OUT=$ROOT/$DST/slurm
mkdir -p "$OUT"

cp "$ROOT/$SRC/conditions.csv" "$ROOT/$DST/conditions.csv"
N=$(( $(wc -l < "$ROOT/$DST/conditions.csv") - 1 ))
PER=$(( (N + W - 1) / W ))
echo "host=$(hostname) src=$SRC dst=$DST alphaP=$ALPHA rows=$N workers=$W per=$PER method=$METHOD"

# PERF=1 samples per-process RSS and CPU time into <campaign>/perf/ (perf_sample.sh).
if [ "${PERF:-0}" = 1 ]; then
  mkdir -p "$ROOT/$DST/perf"
  bash "$S/perf_sample.sh" "$ROOT/$DST/perf/perf_samples.csv" 2 &
  SAMP=$!
fi

PIDS=()
for i in $(seq 0 $((W-1))); do
  START=$(( i * PER ))
  [ "$START" -ge "$N" ] && break
  $PY "$S/sweep_worker.py" --conditions "$ROOT/$DST/conditions.csv" --name "$DST" \
      --start "$START" --count "$PER" --dc-alpha-p "$ALPHA" "${FLAGS[@]}" \
      > "$OUT/w${i}.log" 2>&1 &
  PIDS+=($!)
done
wait "${PIDS[@]}"
[ "${PERF:-0}" = 1 ] && kill "$SAMP" 2>/dev/null
echo "all workers done"
