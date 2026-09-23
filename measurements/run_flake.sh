#!/bin/zsh
# usage: run_flake.sh <treedir> <nruns> <outfile> [testpath] [pytest args...]
#
# Runs a test file N times against a COPY of the tree and records pass/fail per run,
# so `src/` in your worktree can be edited while a measurement is in flight without
# moving the denominator. Make the copy with:
#
#     git archive HEAD | tar -x -C <treedir>
#
# Pass/fail comes from pytest's EXIT CODE, not from parsing its progress line:
#   0 = all selected tests passed (skips are fine)   -> PASS
#   1 = at least one test failed                     -> FAIL
#   other = collection/internal error                -> ERROR, which is not a flake
#
# An earlier version grepped for a line of dots, which silently recorded every run as
# FAIL once `slow` markers started emitting `s` in that line. Do not reintroduce that:
# the progress line is presentation, the exit code is the result.
#
# Two rules for the numbers this produces:
#   - one warm numba cache per tree, throughout; never clear it mid-measurement
#   - the predicate is "one run of the file with at least one failing test", i.e. a
#     per-file-run rate, not a per-test one. Quote it that way.
set -u
TREE=$1; N=$2; OUT=$3
TESTPATH=${4:-tests/test_example_ep_jerk.py}
shift 4 2>/dev/null || shift 3
PY=/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python
cd $TREE || exit 1
: > $OUT
for i in $(seq 1 $N); do
  res=$(PYTHONPATH=$TREE/src $PY -m pytest $TESTPATH -q -p no:warnings --tb=no "$@" 2>&1)
  code=$?
  case $code in
    0) echo "$i PASS" >> $OUT ;;
    1) echo "$i FAIL $(echo "$res" | grep -E '^FAILED' | tr '\n' ' ')" >> $OUT ;;
    *) echo "$i ERROR exit=$code $(echo "$res" | tail -2 | tr '\n' ' ')" >> $OUT ;;
  esac
done
echo "DONE" >> $OUT
