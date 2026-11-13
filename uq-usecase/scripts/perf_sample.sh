#!/bin/bash
# Sample every python and julia process of this user: pid, parent pid, RSS and
# cumulative CPU seconds. Started by qlim_campaign.sh when PERF=1, and killed
# when the workers finish. Summarise with perf_summarise.py.
# usage: perf_sample.sh <out.csv> [interval_s]
OUT=$1
IV=${2:-2}
echo "t,pid,ppid,comm,rss_kb,cpu_s" > "$OUT"
while true; do
  T=$(date +%s)
  ps -u "$USER" -o pid=,ppid=,comm=,rss=,cputimes= \
    | awk -v t="$T" '$3 ~ /^(python|julia)/ {print t","$1","$2","$3","$4","$5}' >> "$OUT"
  sleep "$IV"
done
