#!/bin/zsh
# Re-runs the example do-files of the manual with Stata 17 (-e) and rebuilds the listings in manual/ex/out.
# Usage: NWCORE=/path/to/nwcommands manual/run_examples.sh [ex01_first ex02_effects ...]
# NWDYNAM defaults to the root of this repository; NWPLUS to a fresh scratch directory (so that no installed copy shadows the repository build).
# The data are loaded from the web (nwwebuse) and from manual/ex/data.
here="$(cd "$(dirname "$0")" && pwd)"
export NWDYNAM=${NWDYNAM:-$(dirname "$here")}
export NWPLUS=${NWPLUS:-$(mktemp -d)}
if [[ -z $NWCORE || ! -f $NWCORE/unw_core.do ]]; then echo "set NWCORE to the nwcommands (core) directory"; exit 1; fi
STATA=${STATA:-/Applications/Stata 17.0/StataBE.app/Contents/MacOS/StataBE}
cd "$here/ex" || exit 1
files=("$@"); [[ ${#files} -eq 0 ]] && files=(ex[0-9]*.do) || files=(${^files}.do)
for f in $files; do
  echo "running $f"
  perl -e 'alarm 3000; exec @ARGV' "$STATA" -e do $f > /dev/null 2>&1
done
python3 mklistings.py
