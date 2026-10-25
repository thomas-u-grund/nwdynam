# Gotchas of nwdynam

Known, already diagnosed pitfalls. Entries copied from the nwcommands notes at the time of the split; some mention the core package.

## Mata libraries must be built with the OLDEST supported Stata (Stata 17), never a newer one

The 2026-10-25 releases (core 3.5.0, nwsaom 1.0.7, nwergm 1.2.0, nwdynam 1.1.0, nwrem 1.0.5) shipped `.mlib` files built with Stata 18.
A Stata 17 user got `nwwebuse florentine` -> `<istmt>: 3499 nws_create() not found`, even from a clean install: Stata 17 lists the
library (`mata mlib query`, `mata describe using lnwcommands` work) but silently ignores its contents. The only trace is
`mata which nws_create()`: "(nws_create() in lnwcommands, compiled by Stata 18.0, is too new to be run by this version of Stata and so
was ignored)". The reverse is fine: Stata 18 reads libraries built by Stata 17. So build with the oldest supported Stata (currently 17).
`lib/build.do` of every package stops when `c(stata_version)` is higher than 17 (`NW_ALLOW_NEWER_STATA=1` overrides, for libraries that
are never released). The only reliable check of a release candidate is to run a
command in Stata 17 against the built library. The core has `_nwlibcheck.ado` (called by `unw_defs` and the commands that do not call it) which turns the
failure into a readable message with exit 3499. An OLD library in the user's PLUS folder (for example from a Sep 2026 install) shadows
the new one and gives r(3021) or other class errors; the message tells the user to `ado uninstall nwcommands` first.

## Tied event times are simultaneous (goldfish), not sequential

Before 1.1.0 events with the same time were processed one after another in input order (after the `ordfix` commit: order by time, then row), which gave -652.83 / inertia 5.00 on `socialevolution` with `ceil(t/20)`
while goldfish gives -1595.33 / 3.54. goldfish evaluates all events of a time point at the state before that point (permuting tied rows
does not change its fit; checked 2026-10-25). `DynamState.gfirst` holds the first row of each group and `dynam_group_end()` returns it at the
last event of a group; every likelihood loop updates its tie/lastcontact/tiecount matrices only there. A new loop must do the same, or ties
silently go back to sequential. The native kernel does it with plugin codes 11/12 (v4 = 1 at the last event of a group); the Linux/Windows
plugins are stale until CI rebuilds them. Option `sequential` (scalar `__nwdynam_sequential` read in `DynamState::init`) restores the old behaviour.
nwrem (relevent) still reads ties sequentially. A test that calls a Mata evaluator directly with a scalar left over from a crashed `nwdynam` call can inherit `sequential`; the ado drops it at the end.

## Attaching a covariate to an `eventtime()`-declared network needs `_nwnode`, not a fresh merge

After `nwset sender receiver, eventtime(t) name(netname)` declares an event network (found while
building the `socialevolution` example dataset, 2026-09-12), the dataset left in memory is **not**
the original sender/receiver/t rows anymore — it's one row per actor who appears in at least one
event, with columns `_nwnode` (the actor's own label) and `_nwinclude`. A `merge 1:1 label using
covariates.dta` fails outright (`variable label not found`) because there is no `label` variable to
merge on at that point; the actor identifier is `_nwnode`.

The fix is a rename-merge-rename:

```
rename _nwnode label
merge 1:1 label using covariates.dta, keep(master match) nogen
rename label _nwnode
```

This is unrelated to how covariates work on an already-`nwsave()`d panel dataset (`glasgow`, `s50`,
etc.), where attributes like `smoke1`/`alcohol1` are baked into the saved network file itself and
just show up as ordinary variables on load — no merge needed there at all. The rename-merge-rename
dance is specific to attaching a covariate *after the fact* onto a freshly `nwset`-declared event
network built from a plain, not-pre-declared dataset.

## `lib/build.do` and `deploy_nw_build.do` hardcode the checkout path

Both `cd` into `/Users/tgrund/FILES_NEW/SOFTWARE/nwcommands` by absolute path. Run from a git
worktree or any other copy, they silently rebuild the *main checkout's* mlib instead of the one
you are working on. From a worktree, build with a copy of `lib/build.do` whose paths point at the
worktree.

## A stale installed plugin in PLUS shadows the repo's freshly built native plugin

`findfile()` finds `saom_sim_macos.plugin` (and the other native plugins) in the user's
`~/Library/Application Support/Stata/ado/plus/s/` ahead of `lib/plugins/<platform>/`. After
`make macos-saom_sim` a Stata test can keep running the OLD binary. Don't overwrite the user's
installed copy; copy the new plugin into a scratch dir and `adopath ++ <that dir>` at the top of the
test do-file. `lib/plugins/unix` and `lib/plugins/windows` are built by CI and are stale locally
until it runs.

## A Mata `external` declaration CREATES the external

`external real scalar x` at the top of a function creates `x` (missing) if it does not exist, so a following `findexternal("x") == NULL` test is never true. To test for an optional
external use `pointer(real scalar) scalar p; p = findexternal("x")` without a declaration. (2026-10-23: `SaomBatchSlotTag()` declared `__saom_bslot`, returned "." and every
batch setup sent `BATCHSETUP.|...`, which the plugin sent to its legacy parser: the whole co-evolution suite hung at 100% CPU.) A struct used in a function must be defined EARLIER in the file.

## Copying a rebuilt plugin over a loaded one kills Stata with no message (macOS)

Symptom: `cp new.plugin scratch/lib/plugins/macos/saom_sim_macos.plugin` over an existing copy, then the next
Stata run simply ends in the middle of a command - no error, no `end of do-file` in the log. The crash report
in `~/Library/Logs/DiagnosticReports/StataBE-*.ips` says `SIGKILL (Code Signature Invalid)`.
Cause: macOS caches the code signature of the file; overwriting it in place leaves a stale signature.
Fix: `rm` the old file first, then `cp` (or `cp` to a new name and `mv`). Applies to every scratch directory used for testing.
This is a deliberate
divergence, decided 2026-10-01: the RSiena default looks like a bug (introduced in 1.6.1, commit
1155e31, which moved only outPop/inAct to default 1 and left outPop.c, outPopSqrt, outPopMore,
outPopSqrtMore, outPopThreshold, inAct.c, inActSqrt at 0) and was reported as
https://github.com/stocnet/rsiena/issues/151. Do NOT "fix" nwsaom to match RSiena's default; when
validating against RSiena, set parameter = 1 for these effects. Revisit if the issue is resolved.

## Mata functions and programs defined inside an .ado file are private to that file

`mata: nwergm_estat_bridge_from_netobj(...)` called from a *different* .ado file fails with "not found" (3499): functions
(and `program`s) defined in an .ado are namespaced `file.name` and only callable from programs of the same file. Anything two
helper files need (observed-network bridge, summary functions) belongs in `unw_ergm.do` (the Mata library), as
`ergm_bridge_obs()`. Likewise `estat gof` plots call `_nwergm_gofplot`, a file of its own.

## Stata `syntax` flag options that start with "no" lose the prefix

`syntax [, NOMPLECHECK]` creates the local `mplecheck`, not `nomplecheck`; testing `"`nomplecheck'" == ""` is then always true.
Use a string option (`mplecheck(on|off)`) or a name that does not start with "no".

## `mata mlib index` does not persist: a library installed by net install can be missing from a session's Mata library list

Stata builds its list of Mata libraries (.mlib) once, early in a session; `mata mlib index` rebuilds it for the running session only
and nothing is saved. A library that `net install` put into PLUS/l/ can be absent from the list (`mata mlib query` shows `lnwcommands`
but not `lnwdynam`). The first Mata call then stops with r(3499) "<function>() not found", and a Mata block that mentions a class of
the library stops the whole ado file with r(3000) "class ... undefined (error occurred while loading ...)". Running `mata mlib index`
by hand "fixes" it for that session only, which hid the bug in every manual test.
- Every entry program starts with `_nwdynam_libcheck` (which runs `mata mlib index` once when a library function is not visible); do not
  call Mata before it.
- A **file-scope** `mata:` block that mentions a class of a library would be compiled when the file is auto-loaded, before any program
  runs (see nwergm); neither package has one today.
- Reproduce: `SMOKE_STALE=1 tools/smoke_fresh_install.sh <core> <models>` in the core repo builds exactly this state (the session's Mata
  index is made while only the core library is visible). Without `SMOKE_STALE` a plain new session usually finds all libraries,
  which is why a simple "install, restart, run" test passes.

## Events with the same time are sorted with the row number as second key (fixed 2026-10-25)

`DynamState::init()` sorted the events by time with `order(time, 1)`, which breaks ties by the random-number state: with tied event times the tie matrix evolved in a different order and the log likelihood changed from run to run (socialevolution with time/20 rounded up: -644 to -683 over six seeds). The key is now (time, row number of the network's event list), the order goldfish uses when the times are made distinct in that order. Test: `cscripts/test_nwdynam_ordfix.do`. Note: goldfish treats events with an identical time as simultaneous (its log likelihood on the tied data is -1595.3, not -652.8), which nwdynam does not do; the two agree only for distinct times.
