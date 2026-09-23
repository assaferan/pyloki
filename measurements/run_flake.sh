#!/bin/zsh
# usage: run_flake.sh <treedir> <nruns> <outfile>
TREE=$1; N=$2; OUT=$3
PY=/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python
cd $TREE || exit 1
: > $OUT
for i in $(seq 1 $N); do
  res=$(PYTHONPATH=$TREE/src $PY -m pytest tests/test_example_ep_jerk.py -q -p no:warnings --tb=no 2>&1 | tail -3)
  if echo "$res" | grep -qE '^[.]+\s+\[100%\]'; then
    echo "$i PASS" >> $OUT
  else
    echo "$i FAIL $(echo "$res" | tr '\n' ' | ')" >> $OUT
  fi
done
echo "DONE" >> $OUT
