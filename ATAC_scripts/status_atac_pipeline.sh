#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/lxk/project/lishan-20260613
pid_file="analysis/ATAC/logs/pipeline.pid"
latest_log=$(ls -t analysis/ATAC/logs/run_bulk_atac.*.log 2>/dev/null | head -1 || true)

echo "== ATAC pipeline status =="
if [ -s "$pid_file" ]; then
  pid=$(cat "$pid_file")
  ps -p "$pid" -o pid,ppid,stat,etime,cmd || true
  echo "-- child processes --"
  pgrep -P "$pid" -af || true
else
  echo "No PID file found."
fi

echo "-- outputs --"
for d in trimmed bam peaks bigwig counts annotation figures matrix qc logs; do
  if [ -d "analysis/ATAC/$d" ]; then
    printf "%-12s %s files\n" "$d" "$(find "analysis/ATAC/$d" -type f | wc -l)"
  fi
done

if [ -n "$latest_log" ]; then
  echo "-- latest log: $latest_log --"
  tail -80 "$latest_log"
fi
