#!/bin/bash
set -e
THRESHOLD_PAGES=1000000  # ~16GB free (16384-byte pages)
echo "Waiting for memory headroom (need >${THRESHOLD_PAGES} free pages)..."
while true; do
  free_pages=$(vm_stat | awk '/Pages free/ {gsub("\\.", "", $3); print $3}')
  echo "$(date '+%H:%M:%S') free_pages=${free_pages}"
  if [ "${free_pages}" -gt "${THRESHOLD_PAGES}" ]; then
    echo "Memory headroom OK, proceeding."
    break
  fi
  sleep 60
done

source /Users/assaferan/Documents/GitHub/pyloki/.venv/bin/activate
cd /Users/assaferan/Documents/GitHub/pyloki
rm -rf /Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment/results_economized
python /Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment/run_search.py economized
echo "WAIT_AND_RUN_DONE"
