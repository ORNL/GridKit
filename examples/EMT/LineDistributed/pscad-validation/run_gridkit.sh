#!/usr/bin/env bash
# Run every GridKit suite in scenarios/gridkit with the line_distributed executable.
#
#   run_gridkit.sh <line_distributed executable> [suite name pattern ...]
#
# Results go to results/gridkit/<line>__<model>/<scenario>/ (input.json, samples.csv,
# steps.csv, statistics.json), as the executable writes them. Intended for the
# gridkit-sdk container, for example:
#
#   docker exec gridkit-ld bash -lc '/work/examples/EMT/LineDistributed/pscad-validation/run_gridkit.sh \
#       /home/gridkit/build/examples/EMT/LineDistributed/line_distributed'
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]:-.}")" && pwd)"
executable="${1:?path to line_distributed}"
shift || true
mkdir -p "$here/results/gridkit"
summary="$here/results/gridkit/runs.json"
echo "{" > "$summary.tmp"
first=1
for suite in "$here"/scenarios/gridkit/*.json; do
  name="$(basename "$suite" .json)"
  if [ "$#" -gt 0 ]; then
    match=0
    for pattern in "$@"; do [[ "$name" == *"$pattern"* ]] && match=1; done
    [ "$match" -eq 1 ] || continue
  fi
  echo "== $name"
  start=$(date +%s.%N)
  status=0
  # A pathological run (see README: undamped recirculating fronts) must not block the rest.
  timeout 1800 "$executable" "$suite" "$here/results/gridkit/$name" > "$here/results/gridkit/$name.log" 2>&1 || status=$?
  cat "$here/results/gridkit/$name.log"
  [ "$status" -eq 0 ] || echo "FAILED $name (exit $status)"
  end=$(date +%s.%N)
  [ "$first" -eq 1 ] || echo "," >> "$summary.tmp"
  first=0
  printf '  "%s": {"seconds": %s, "exit": %s, "finished": "%s"}' "$name" "$(awk "BEGIN {printf \"%.3f\", $end - $start}")" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$summary.tmp"
done
echo "" >> "$summary.tmp"
echo "}" >> "$summary.tmp"
mv "$summary.tmp" "$summary"
echo "DONE"
