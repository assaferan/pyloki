#!/bin/bash
# Memory-aware launcher: wait until enough RAM is free, then run the search.
#
# Written when magma.exe processes from unrelated sessions were consuming 6-18 GB
# each and OOM-killing runs. Check `ps aux | grep magma` / `vm_stat` before
# assuming a failure is pyloki's fault.
#
# Paths are derived from this script's location so it works in any worktree.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PY="$REPO/.venv/bin/python"

NEED_GB=16
while true; do
  FREE_GB=$(vm_stat | awk '/page size/{ps=$8} /Pages free/{f=$3} /Pages inactive/{i=$3} END{printf "%.0f", (f+i)*ps/1073741824}')
  if [ "$FREE_GB" -ge "$NEED_GB" ]; then
    echo "$(date '+%H:%M:%S') ${FREE_GB}GB free >= ${NEED_GB}GB, launching"
    break
  fi
  echo "$(date '+%H:%M:%S') only ${FREE_GB}GB free, waiting 60s"
  sleep 60
done

rm -rf "$HERE/results_economized"
"$PY" "$HERE/run_search.py" economized
