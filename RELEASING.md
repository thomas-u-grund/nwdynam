# Releasing nwdynam

## Master procedure: https://github.com/thomas-u-grund/nwcommands/blob/main/RELEASING.md

All release procedures and policies of the five packages (steps, order, test launch, tags, GitHub release, website update, Zenodo, lockstep with the core
Mata library, building libraries with Stata 17, known traps) are in the master document above. Follow it; this file only holds what is specific to nwdynam.
Two rules from there that are easy to miss: push the release tag by name (`git push origin vX.Y.Z`), never `git push origin main --tags`
(old local tags of this clone never existed on GitHub); and build `lib/*.mlib` with Stata 17, not 18.

## Version files and checks

- Single source of truth: `VERSION` (line 1 the version, line 2 the release date in Stata format, e.g. `25oct2026`). Semantic versioning (master B.4, C4).
- `python3 tools/set_version.py X.Y.Z [date]` writes it into every display (first line of every `.ado`, help headers and the footer of every help file, `nwdynam.pkg`, `stata.toc`,
  `README.md`, `NOTICE.md`, `CITATION.cff`, the `e(version)` / `nwdynam version` literals); `--check` changes nothing.
- Add `## [X.Y.Z] - YYYY-MM-DD` to `CHANGELOG.md` (Added / Changed / Fixed / Removed); its text is the GitHub release notes.
- `NWCORE=/path/to/nwcommands python3 tools/check_version.py` = version check plus the core-library hash check (below).
- Install lines: latest `net install nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main/) replace`; a given release `.../nwdynam/vX.Y.Z/`.
  `nwdynam.pkg` states `Requires nwcommands (core), version A.B.C or later`; raise it (and `README.md`) only when this release needs a newer core.

## Library and lockstep with the core

`lib/build.do` builds `lib/lnwdynam.mlib` from `unw_dynam.do`; it does not load the core and the sources do not mention the core class `nw_def` (verified 2026-10-25), so a changed core library cannot cause the r(3021) class mismatch here.
The package nevertheless carries `lib/core_mlib.sha256`, and `NWCORE=... python3 tools/check_version.py` fails when the core library's hash differs. In a lockstep batch (master C3) the practice so far is: rebuild with Stata 17 (`<Stata 17> -e do lib/build.do` from the repository root), run
`NWCORE=/path/to/nwcommands python3 tools/record_core_lib.py`, test and release a patch. Whether a record-only update without a release is enough is an open decision (master B.5).
- `nwinstall, buildmata(nwdynam)` (core `nwinstall_buildmata.ado`) rebuilds lnwdynam.mlib on the user's computer from `unw_dynam.do` on `main`. A new or renamed Mata source file or library of this package must be added to `nwinstall_buildmata_tab` in the core (see its GOTCHA.md), and the `mata mlib add` patterns there must match `lib/build.do`.

## Tests

`NWCORE=/path/to/nwcommands cscripts/run_tests.sh [jobs]` (first step: `tools/check_version.py`; each test runs in its own `-e` Stata process; set `STATA` to the Stata 17 binary, then repeat with 18).
Most tests print no PASS marker: failure is an `r()` error and a missing log is NOT RUN (master C0).

## Plugin

`native/dynam_sim.c` builds `dynam_sim` (CI commits the binaries to `lib/plugins/` on `main`; `pull --rebase`, tag after the bot commit; nwdynam v1.1.0 is tagged on the bot commit). Ties in event times are simultaneous as in goldfish (`sequential` for the old way); events with equal time are sorted with the row number as second key.
A stale installed plugin shadows the repo build: scratch `adopath`, `rm` before `cp`.
Every fit needs the plugin (no Mata fallback). The protocol number is `DYNAM_PROTOCOL` in `native/dynam_sim.c` and `DYNAM_PROTOCOL_NEEDED()` in `unw_dynam.do`; raise both together whenever the blob layout or the meaning of an argument changes, so that an old plugin is recognised (r(3495)) instead of giving wrong numbers. A release rebuilds `lib/lnwdynam.mlib` with Stata 17 and the three plugins (macOS here, Linux and Windows by CI) from the same commit: until CI has rebuilt Linux and Windows those platforms stop with r(3495). The Mata evaluators of the earlier versions are test-only reference code (`cscripts/ref_dynam_mata.do`); `test_nwdynam_plugin_{rate,choice,coord}.do` compare the plugin with them. Protocol 3 (nwdynam 1.6.0, time-varying covariates) added the change-event block to the blob; the Mata references do not know it, so `test_nwdynam_timevarying.do` compares the plugin with goldfish instead (data and reference values in `dev/tvc/`, made by `dev/dynam_tvc_crosscheck.R`; regenerate with `Rscript dev/dynam_tvc_crosscheck.R dev/tvc`, goldfish 1.6.12). The Linux and Windows plugins of 1.5.1 report protocol 2 and stop with r(3495) until CI has rebuilt them; check that after the CI commit.

## Manual

There is no nwdynam manual and no Zenodo record; the package is cited by repository and version. No website "Manual" menu entry until one exists.
