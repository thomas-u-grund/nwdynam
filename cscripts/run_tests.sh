#!/bin/zsh
# nwdynam test runner (macOS). Usage: NWCORE=/path/to/nwcommands cscripts/run_tests.sh [jobs]
# NWCORE: directory of the core nwcommands package (unw_core.do and the core .ado files).
# The tests are independent do-files; each runs in its own Stata process (-e); logs are run_<test>.log in the repo root.
cd "$(dirname "$0")/.." || exit 1
python3 tools/check_version.py || exit 1
if [[ -z $NWCORE || ! -f $NWCORE/unw_core.do ]]; then echo "set NWCORE to the nwcommands (core) directory"; exit 1; fi
export NWCORE
jobs=${1:-3}
STATA=${STATA:-/Applications/Stata/StataBE.app/Contents/MacOS/StataBE}
tests=(test_nwdynam_ado test_nwdynam_native test_nwdynam_twomode_ado test_nwdynam_fixes1 test_nwdynam_ordfix test_nwdynam_ties)
tmp=$(mktemp -d)
run_one() {
  local t=$1
  echo "do cscripts/$t.do" > $tmp/run_$t.do
  $STATA -e do $tmp/run_$t.do > /dev/null 2>&1
  printf '%-34s PASS lines: %3s   errors: %s\n' $t "$(grep -c -E 'PASS|ALL OK|ALL TESTS PASSED' run_$t.log)" "$(grep -E '^r\([0-9]+\)' run_$t.log | sort -u | tr '\n' ' ')" > $tmp/res_$t.txt
  rm -f run_$t.log
}
running=0
for t in $tests; do
  run_one $t &
  running=$((running + 1))
  if (( running >= jobs )); then wait -n 2>/dev/null || wait; running=$((running - 1)); fi
done
wait
for t in $tests; do [[ -f $tmp/res_$t.txt ]] && cat $tmp/res_$t.txt; done
rm -rf $tmp
echo "done"
