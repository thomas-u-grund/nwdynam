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

## A `/*` inside a `*` comment line opens a block comment (found 2026-10-26)

A do-file comment such as `* a tree without lib/plugins/* (...)` contains `/*`, which starts a block comment that swallows the rest of the file: the log shows the whole file echoed and nothing runs, with no error. Never write `/*` inside a `*` comment line.

## `nwinstall, buildmata(nwdynam)` compiles the source files listed in the core (2026-11-03)

The user-side rebuild of lnwdynam.mlib lives in the core (`nwinstall_buildmata.ado`, table `nwinstall_buildmata_tab`): it lists the source files of this package, the library names and the `mata mlib add` patterns copied from `lib/build.do`. If this package gets a new Mata source file (as nwsaom did with `unw_saom_struct.do`, included by `unw_saom.do`) or another library, the table must change, or the rebuild fails with "file ... not found" (r(601)) or leaves functions out. The sources fetched are those on `main`.

## Mata `_error("long text")` loses the message (r(3300))

`_error(errtxt)` fails itself with r(3300) "argument out of range" when the text is longer than about 110 characters,
and the intended message is never shown. Print long messages with `errprintf()` and leave with `exit(rc)`.

## The native engines switch to the frame `__nwdynam_native`

`DynamFitUnit1()`/`DynamFitRateUnit1()` work in a frame of their own and switch back only on the normal path. An error
inside Mata (for example the optimizer failing) used to leave the user in that frame, and every later call failed with
r(3598) "frame already defined". `nwdynam.ado` therefore runs the fit inside `capture noisily { }` and restores the
frame, drops the helper frame, scalars and Mata globals before passing the error on.

## Event lists can hold actor index 0 or a missing time

The core's `build_eventlist()` gives index 0 to an event whose sender or receiver label is not a node name (missing id,
or numeric ids >= 1e7, labelled with %9.0g in core nwset), and `nwset` keeps missing event times. The native engines
do not check either; `nwdynam.ado` refuses both before the fit (r(416)).

## Replay without arguments

`nwdynam` alone used to fit the default model of the current network, and `estimates replay` ran exactly that, so it
refitted. Now `nwdynam` alone replays when `e(cmd)` is nwdynam, and still fits otherwise.

## Event values: column 4 of the event list, the native kernel only gets them when some value is 0 or below

Since the core keeps event values (`nwset sender receiver value, eventtime(t)`), `get_eventlist()` returns 4 columns
(sender, receiver, time, value; value 1 without a value variable); `DynamState::init` pads an older 3-column list with 1.
Every likelihood loop keeps `tieval` (summed values, goldfish's increments): unweighted effects test `tieval > 0`
(inertia, recip, closure effects), but the unweighted indeg/outdeg add `sign(tieval)` as goldfish does (a tie below 0
counts -1; inside a window this happens when a negative event outlives the event it cancelled - found because the rate
window test disagreed with goldfish while the choice one agreed). Windows keep their own summed matrix, expired with
`dynam_window_expire()`. With only positive values the binary state equals the old "any event" state, so the native
kernel gets the values (v5, codes 101/102/111/112) only when some value is 0 or below: the stale Linux/Windows plugins
answer those codes with r(198) until CI rebuilds them.

## optimize() cannot polish a Newton fit to round-off (rate intercept, 1.4.0)

`optimize()` accepts a step only if the log likelihood rises, and at a log likelihood of -100 (about 1e-14 resolution) the last 1e-8 of the parameters is below that resolution; at -400000 it is far worse and the iterations end in endless "(backed up)" lines (a 50,000-event fit then ran for minutes through the BFGS and old-engine fallbacks). `dynam_rateint_newton()` therefore lets `optimize()` stop at a Newton decrement of 1e-11 and finishes with explicit Newton steps judged by the gradient norm. Do not tighten `conv_nrtol`/`ptol`/`vtol` instead: that only produces more backed-up iterations. The result is accepted through `dynam_rateint_accept()` (decrement below 1e-8), never through `optimize_result_converged()`.

Mata `exp()` of a number above about 709.78 is missing and `sum()` ignores missing values: an evaluator that sums `exp(lp)` silently loses the hazard term and a line search then sees a spurious improvement. Return a missing log likelihood instead (the cached and the native evaluator do, above lp = 700).

The plugin code 3 (rate with intercept, indeg/outdeg) is absent from plugin binaries of earlier releases: the fit probes the loaded plugin once with `_stata(..., 1)` and falls back to the cached Mata evaluator when the call fails (tested with the 1.3.2 macOS plugin). The Linux and Windows plugins are rebuilt by the CI workflow on a push that touches `native/`.


## Plugin protocol: an old plugin crashes Stata on a short argument string, and a loaded plugin cannot be replaced

nwdynam 1.4.0 and earlier plugins (protocol 1) read three numbers from the argument string BEFORE they look at the algorithm code; a probe such as
`"99"` makes `atof(NULL)` crash Stata (no error message, the log just stops). `DynamNat::open()` therefore probes with `"99 0 0 0 0 0 0"`; the old plugin answers r(198) and the glue
reports r(3495). Within one Stata session the first plugin loaded for the program name `dynamnative` stays, even when another file of the same name is put first on the adopath
(that is why `test_nwdynam_oldplugin.do` and `test_nwdynam_noplugin.do` are separate files). Stata's `shell` passes compound quotes (`` `"..."' ``) through to the shell: use plain `"..."` in shell lines. `git -C <dir>` works in `shell`, `cd dir && ...` did not.

## The Mata evaluators are reference code only

`cscripts/ref_dynam_mata.do` holds the former evaluators (`dynam_*_loglik_grad_*`, the cached rate-with-intercept evaluator); `unw_dynam.do` has none. A test that needs them runs `do unw_dynam.do` and then `do cscripts/ref_dynam_mata.do`.
The plugin tests compare log likelihood, gradient and Hessian on random cases (`test_nwdynam_plugin_*.do`); when a statistic changes, change it in BOTH places. The fit tests use `lib/lnwdynam.mlib`: rebuild it with Stata 17 (`lib/build.do`) after every change of `unw_dynam.do`, or the test silently runs the old Mata.

## optimize()'s "nr" technique backs up for hundreds of evaluations at round-off level (found 2026-11-16)
On dense data (large log likelihood) the Newton decrement stops falling at about 1e-11 and Mata's `optimize()` with technique `nr` keeps taking "backed up" steps (up to 40 iterations times 30 to 50 evaluations each) because it never meets its tolerances. The estimates were right, the time was 20 to 100 times too long (`dynam_nat_optimize` mode 0/2 and the BFGS polish). `dynam_rateint_newton()` therefore runs its own damped Newton loop on the exact Hessian and stops at a decrement below `max(1e-12, 1e-13*|ll|)`. Do not go back to `optimize_init_technique("nr")` for these fits; count evaluations (a global counter in `dynam_nat_eval`) when a fit is slow before suspecting the kernel.

## Time-varying covariates: what goldfish does, so that the plugin does the same

(Plugin protocol 3, nwdynam 1.6.0.) The facts below were found by comparing the plugin with goldfish 1.6.12 on random logs
(`dev/dynam_tvc_crosscheck.R`, test `cscripts/test_nwdynam_timevarying.do`); they are easy to get wrong:
- A change event at time t is applied AFTER the dependent events of time t (goldfish takes the dependent list first at equal times): the plugin applies the changes with time < t before it evaluates a group of events at time t.
- Rate with `intercept`: goldfish's observation period starts at the earliest event of ANY list and ends at the latest. A change before the first dependent event adds the exposure from that change to the first event (the first event then has a positive waiting time), a change after the last event adds the exposure up to it. Every change inside a waiting time splits it (`eval_rateint`). Choice, rate without intercept and coordination are not affected by changes outside the log.
- Several changes of one actor at one time: the later row wins (increments add up in row order, done in `dynam_change_rows()`).
- The change network is an ordinary `nwset` event network (sender = receiver); its value variable is the NEW value (the dependent network's value is an increment). Without a value variable nwset stores 1 for every event; the plugin cannot tell, so every change would set 1.
- The plugin gets the changes as absolute new values (the Mata glue turns increments into values) in blob block "change events"; tertius keeps its running numerator (`tnum`), which a change of the covariate updates for the actor's out-neighbours.
- Coordination standard errors come from a numerical derivative of the gradient (no analytic Hessian) and differ from goldfish's by up to a few percent, with or without changes; the cross-check allows 10% for them.
- The old Mata evaluators (`cscripts/ref_dynam_mata.do`) do NOT know time-varying covariates (owner's decision of 2026-11-16: validate against goldfish instead).

## Composition change: what goldfish does, so that the plugin does the same

(Plugin protocol 4, nwdynam 1.7.0.) Found by reading goldfish's source (`DyNAM_*_default.cpp`, `convert_rdata_to_cdata.cpp`, `class_checks.R`) and
comparing with goldfish 1.6.12 on random logs (`dev/dynam_comp_crosscheck.R`, test `cscripts/test_nwdynam_composition.do`):
- Presence changes with time <= t are applied BEFORE the event at time t (`C_convert_composition_change` counts the changes with time <= the reference time),
  the opposite of covariate changes (those with time < t). So a join at the time of an event lets the actor take part, a leave at that time excludes it.
- Choice: the candidate receivers are the present actors (not the sender). Rate: the hazard / normaliser sums over the present senders only. Coordination: pairs
  with an absent member get probability 0 (the row normalisation uses the present partners; an absent actor's own row is all zero). Statistics are computed
  for everybody (no effect function looks at `present`).
- Rate with intercept: the composition changes are NOT events of their own: they do not split waiting times. The waiting time that ends at a reference
  row (a dependent event, or a covariate change when `change()` is used) is exposed to the actors present at that row's time. A leave inside the waiting
  time removes the actor from the whole waiting time, a join at its end adds it to all of it. A composition change before the first event does not start
  the observation (only covariate changes do).
- An event of an absent sender or receiver is an error in goldfish (`checkPresence`, "the node is not present at time"); so is a covariate change event of an
  absent actor (when the composition list is linked first). nwdynam checks both in `dynam_check_presence()` (Mata, before the fit; the plugin refuses an absent
  sender/receiver as well). Two-mode: sender (mode 1) and receiver (mode 2) are both checked.
- `present(var)` is only a convenience: it equals a leave event before the first event. opportunities() cannot be combined with present() (data layout), use the leave event.
- Absent actors that are not in the log at all are not actors of an `nwset` event network: add them with `nwaddnodes` (one-mode) before the fit; in a two-mode network
  `nwaddnodes` appends the new node at the END of the node list, which breaks nwdynam's "mode 1 first" assumption (r(198), "sender must be mode 1"): not solved here,
  the two-mode cross-checks only contain actors with events.
- Times that must be equal (a join at the time of an event) must be stored with the same precision: an imported `float` event time 63.31 is 63.3100013..., a `double`
  63.31 is smaller, so the "same time" change would be earlier. Use float variables for such changes or read both from the same source.
- The old Mata evaluators do NOT know composition change (owner's decision: no new Mata references); the goldfish files in `dev/comp` are the reference.
- Windows/Linux plugins in the repository are stale (protocol 3) until CI rebuilds them: those platforms stop with r(3495).

## Windows on tertius, coordination windows, two-mode effects, opportunities with covariates: what goldfish does

(Plugin protocol 5.) Data and reference values: `dev/win/` (made by `dev/dynam_win_crosscheck.R`, regenerate with `Rscript dev/dynam_win_crosscheck.R dev/win`,
about 15 minutes; test `cscripts/test_nwdynam_combinations.do`).
- Event times in the core are stored as float: 125.1 is 125.0999985. A window expires when `time_of_tie + window < current time` (the sum computed in double as in R),
  so a boundary case is decided by the stored float. The cross-check logs use times that are multiples of 0.25 (exact in float) and are read back with `recast double`.
  `import delimited` also makes decimals float: recast to double before `nwset`, or the boundary cases differ from goldfish.
- Windowed rate effects with `intercept`: the hazard integral is split at the expiry of every window as well as at the covariate changes. Expiries after the last
  change/event do not count. A creation at time t is applied before a dissolve at the same time.
- Windowed tertius (goldfish `update_REM_choice_tertius`): a per-receiver cache holds the mean attribute of the in-neighbours in the windowed network (NaN if none). When a tie
  flips sign: if all caches are 0, or any is NaN, or the new value is NaN, every NaN-cache receiver gets the mean of the non-NaN caches as its statistic. goldfish turns
  a missing statistic into 0 (found with `preprocessingOnly = TRUE`).
- A change of the tertius variable together with `tertiuswindow()` and `intercept` (rate) is refused: goldfish's likelihood there differs from every ordering of the
  change and the window expiry tried (about 1 log-likelihood unit on some logs; the statistics of the plugin and of goldfish's own preprocessing agree, the difference is in the
  start and end of the observation). Do not try to lift this without a new idea.
- choice_coordination: a windowed and a weighted inertia (or indeg) in one model are nearly collinear and the likelihood is flat: the estimates agree to about 1e-4, not 1e-6,
  although the log likelihood agrees to 1e-6. The standard errors are numerical (10% tolerance). The generator skips such combinations.
- Two-mode (goldfish 1.6.12): `alter(isTwoMode = TRUE)` and `ego(isTwoMode = TRUE)` need the flag; `egoAlterInt` has no flag and takes mode-1 attribute x mode-2 attribute; `tertius` uses
  the mode-1 attribute; `four` works. `same`, `diff`, `sim` stop with an error; `commonSender` returns an n2 x n2 matrix; recip, outdeg (choice), trans, cycle, commonReceiver, nodeTrans
  and (rate) indeg, tertius, nodeTrans are refused by goldfish. Actor labels sort as text: `A10` comes before `A2`, so an `opportunities()` actor index differs from the
  generator index when there are 10 or more actors (the cross-check uses at most 9).
- `opportunities()` takes one row per event and available actor, the covariates one row per actor: `oppframe(name)` reads the two variables from a frame.


## Effects on other networks: what goldfish does (netterms, ignorerep)

- Statistics are literal translations of goldfish's update functions: `d = sign(new) - sign(old)`, s and r skipped, caches start at 0 for windowed terms but keep the start matrix. Self-ties in a static network count as in goldfish (`nwset, mat() selfloop` in tests; without `selfloop` nwset drops the diagonal).
- Order at equal times: dependent events are evaluated first, then exogenous events; creations before expiries; network 1 before network 2.
- `four` recomputes the paths from the network after the change (goldfish `update_DyNAM_choice_four`).
- tertiusDiff: goldfish stops with an error ("subscript out of bounds", "cannot update") when a windowed term leaves a receiver without in-neighbour; nwdynam stops too (only for windowed terms; without a window goldfish copes with removals). goldfish's matrix-wide mean imputation of NaN in `getIterationStepState` is not implemented. Not available: time-varying attributes, two-mode.
- `ignorerep` in the rate model: goldfish accepts it but no hypothesis reproduced its likelihood, so it is refused.
- Refused because goldfish does not offer them: cycle/commonsender/commonreceiver/mixed cycle and common effects in choice_coordination; `weighted` on closures.
- Plugin protocol stays 5: the term block is appended to the blob (blob[40] = its length); the Mata glue passes it in the external `__nwdynam_xblob`.

## Postestimation: the fit frame, goldfish's record count, ties in Top

* A fit keeps its data in the frame `__nwdynam_fit` (the plugin blob in v1, plus `ord` and `nm`; characteristics `bloblen`, `model`, `fitid`) instead of closing it; `e(fitid)` ties the estimates to it, `estat` stops with r(498) otherwise. Fits do no extra work. The next fit replaces the frame; `clear all` removes it.
* `estat` only reaches `nwdynam_estat` when `e(estat_cmd)` is set; without it Stata runs its own `estat ic`.
* The new plugin codes 30 (per-event values) and 40 (outliers, change points) write their output behind the blob in v1 and leave the blob and codes 20/99 unchanged: no protocol bump. A plugin without them stops with r(3495).
* Rate with intercept: goldfish's `nEvents`, `intervalLogL` and `eventProbabilities` also hold the right-censored intervals (`orderEvents == 2`) that covariate changes and window expiries create; AIC and BIC use that count. `examineOutliers` / `examineChangepoints` subset to the events. Composition changes do not split periods. The "probability" returned for this model is `lp_j - timespan * sum(rates)`, not a probability.
* Per-event values are only comparable to goldfish at the same coefficients: tests repost goldfish's `e(b)` before `estat` (the fit agrees only to about 1e-7).
* `examineOutliers` default Hampel window is `N/2 - 1` (the loop `(w+1):(N-w)` then tests only the middle event(s); an even N with 2w+1 > N errors). `method(top)` ranks by `order(x)`; exact ties in the log likelihood (events with identical state) are broken by rounding noise, so goldfish may name another event.
* `changepoint`: minseglen defaults to the largest run of equal times; the variance cost needs minseglen >= 2 and N >= 4; BinSeg errors when Q = 5 > N/minseglen + 1; the penalty is MBIC, `3 log n`. The list of change points always ends with n; goldfish removes entries with equal times (`!duplicated(time[cpt], fromLast = TRUE)`).
* Mata: one-line brace bodies (`if (x) { a; b }`, `} else { ... }` on one line) gave "illegal arglist" in this file; use multi-line blocks. `st_numscalar()` of a missing scalar aborts; test with `_stata("confirm scalar ...", 1)`.

## Mata passes arguments by reference; coordination needs a pre-optimisation (found 2026-11-16, fixed in 1.8.1)

- A Mata function that sorts, re-indexes or edits one of its arguments changes the CALLER's object. `nwdynam.ado` passes `*(netobj->get_eventlist())` and `*(netobj->get_matrix())` of the network behind `netterms()` straight into `dynam_xlayer()`, which used to sort the event list and replace the actor indices in place: the first fit was right, every later fit with the same event network re-mapped the already mapped indices (log likelihoods between -2272 and -2307 instead of -2297.04, a term near -28). Copy such arguments into locals first. The same holds for any new glue that hands a network's own matrix or list to a function.
- choice_coordination (no Hessian, BFGS in Mata): the first BFGS step is the whole gradient. With a statistic of large values (indeg counts up to about 100) it overshoots, exp() overflows, the log likelihood is missing and `optimize()` reports "flat region" (r(111)); the later restarts do not help. `dynam_nat_warm()` runs a damped quasi-Newton ascent (backtracking, missing = failed step) before each start, `optimize()` then polishes and gives the variance. goldfish: indeg alone 0.09896, ll -5250.96 (first 631 emails as undirected contacts).
- `estat changepoints, plot` used the matrix `idx` after `return matrix cpts = idx` had moved it away; compute everything that needs a returned matrix before the `return matrix` line.
