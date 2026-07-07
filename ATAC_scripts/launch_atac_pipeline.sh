#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/lxk/project/lishan-20260613
mkdir -p analysis/ATAC/logs
nohup bash -lc 'THREADS=20 analysis/ATAC/scripts/run_bulk_atac.sh' \
  > analysis/ATAC/logs/launcher.nohup.log 2>&1 &
echo $! > analysis/ATAC/logs/pipeline.pid
echo "started_pid=$(cat analysis/ATAC/logs/pipeline.pid)"
