# Changelog

All notable changes to nwdynam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.11.2] - 2026-11-26

### Fixed
- Help: the examples declare their event networks as `directed`; without it, `nwset` made an undirected network and the examples stopped.
- Tests: `nwset` with `dropvars` where unused variables are discarded (needed to run the tests with core 3.39.0; the package itself works with core 3.36.0 or later).

## [1.11.1] - 2026-11-22

### Changed
- Rebuilt for nwcommands 3.36.0; requires core 3.36.0 or later. No other change.

## [1.11.0] - 2026-11-21

### Changed
- `estat choiceplot`: ties in the fitted probabilities are handled in one stated way. The rank of the realised receiver (pair), and so the mean and the median rank, is the mid rank of its tie group. A top-k share counts, for an event whose realised member ties with others for places g+1 to g+m, the fraction (k-g)/m (between 0 and 1) of the places up to k: the exact expectation of the share under a random order within the tie group. It is deterministic (no random numbers) and differs from counting the mid rank as a hit only for events whose tie group straddles k (manufacturing emails, k = 1: .1480, a mid-rank count would give .1432). The numbers of 1.10.0 are unchanged; the help now says which rule applies to which statistic (the 1.10.0 note "ties share the mid rank" holds for the rank only). goldfish has no such statistic, so there is no reference convention to follow.
- When ties occurred, a line under the table gives the number of events whose realised member ties and the rule; new returned scalar `r(N_ties)` and new variable `nties` (size of the tie group) in the frame of `frame()`.
- The estimates and every other command are unchanged; the compiled plugin is the one of 1.9.1 (the Mata library is rebuilt). New test `cscripts/test_nwdynam_choiceties.do`.

## [1.10.0] - 2026-11-20

### Added
- `estat choiceplot` (choice and choice_coordination): the rank of the realised receiver (pair) among the risk set of every event by the fitted probability. Graph of the share of events with rank at most k against the uniform expectation and the expectation under the fitted model; mean and median rank, and a table of top-k hits (`topk()`, `maxrank()`, `frame()`, `noplot`, `name()`, any `twoway` option). Ties share the mid rank.
- `estat rateplot` (rate, with or without `intercept`): observed against expected number of events per actor as sender with the 45 degree line (`over(actor)`), and observed and expected events per period for the most active actors (`over(time)`, `periods()`, `topactors()`, any `twoway` option; count axes have integer labels). For the model with intercept the expected counts are the shares of the total hazard and add up to the observed number of events.
- Both are extensions; goldfish has no such plots. They use the probabilities of the native evaluator (no Mata likelihood). Test `cscripts/test_nwdynam_resultplots.do` checks the numbers against an independent enumeration. The Mata library gains `dynam_post_rank()` and `dynam_post_rate()`; the plugin is unchanged.
- Graph options of both: default graph names `choiceplot` and `rateplot` (`graph display choiceplot` works without `name()`); `name(name)` accepts a trailing `, replace`; `saving()` and `export()` as above; default titles and axis titles are black and `title("")` suppresses the default title; a mistyped graph name, a bad `saving()`/`export()` file or directory stops with a message before any computation. With `noplot`, `name()`, `saving()`, `export()` and `twoway` options give return code 198.
- The tables and returned results of the existing diagnostics (`estat ic`, `eventlik`, `eventprob`, `outliers`, `changepoints`) are unchanged.

### Changed
- The graph options of the postestimation plots follow the common convention of the nw packages and use the graph helper of the core: **nwcommands (core) 3.29.0 or later is required** (asked when a graph is drawn; with an older core the plot stops with r(199) and a message, tables and all other commands are unaffected).
- `estat outliers, plot` and `estat changepoints, plot`: new options `saving(filename)` (a `.gph` file, replaced if it exists, `, replace` accepted), `export(filename [, graph_export_options])` (`graph export` of this graph) and any `twoway` option, passed last so that it overrides the defaults. Options without `plot` stop with return code 198, except `name()` alone, which is still accepted. The default graph is unchanged (unnamed, the picture is byte-identical to 1.9.1).

Tested on macOS with Stata 17 BE (plot tests, the ado test and the result-plot test, against core 3.29.0); on Linux and on Windows with StataNow 19.5 MP (the plot tests and the package tests of this version). Stata 18 on the Mac was not run for this version. The compiled plugin is the one of 1.9.1 (unchanged).

## [1.9.1] - 2026-11-18

### Fixed
- On a Mac, the console version of StataNow (19.5) reports `c(os)` as `Unix`, so the package looked for the Linux plugin and stopped with an error that the plugin was too old or could not be loaded. The platform is now taken from `c(machine_type)` as well, so a Mac always uses the macOS plugin.

## [1.9.0] - 2026-11-17

### Changed
- Coefficient tables show the coefficient names in full by default: the confidence intervals are left out, as with `noci`, at a `linesize` below 100 (where names of up to 32 characters do not fit next to them); before, long names were cut (e.g. `weightedine~a`). The new option `ci`, or `level(#)` given explicitly, shows the intervals as before. The display options of Stata's estimation commands (`noci`, `nopvalues`, `cformat()`, `pformat()`, `sformat()` and the others of `ereturn display`) are accepted in fits and in the replay (`nwdynam` alone, `estimates replay`). Tables whose names all fit, and all output outside the coefficient tables, are unchanged.
- `ci` together with `noci` stops with return code 198; `noci` is accepted.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (full test suite). Windows and Linux: not tested for this version; the compiled plugin is the one of 1.8.0 (unchanged).

## [1.8.1] - 2026-11-16

### Fixed
- `netterms()` with a network declared with `eventtime()`: the command re-indexed and sorted the network's own stored event list in place. The first fit was correct, but every later fit in the same session with the same event network gave wrong estimates, and sometimes an apparently unidentified term. Fits now work on a copy and are identical on every call. Earlier results are affected only if they came from a second or later fit with an event network whose actor order differs from that of the fitted network.
- A static network in `netterms()` no longer has its missing cells overwritten by 0 in the stored network.
- `estat changepoints, plot` stopped with "not found" and r(100). Fixed.
- `submodel(choice_coordination)` with statistics of large values, for example `indeg`, stopped with "flat region", r(111), although the model is identified. Coordination fits now run a damped quasi-Newton pre-optimisation before the final optimisation; the estimates agree with goldfish 1.6.12.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (full test suite, including a new test for these fixes). Windows and Linux: not tested for this version; the compiled plugin is the one of 1.8.0.

## [1.8.0] - 2026-11-16

### Added
- **More combinations, as in goldfish.** Windowed rate effects with `intercept` (`indegwindow()`, `outdegwindow()`); a window on `tertius` (new option `tertiuswindow(#)`) in the choice, rate, rate-with-intercept and coordination models; `inertiawindow()`, `indegwindow()`, `weightedinertia` and `weightedindeg` in `choice_coordination`; in two-mode networks the effects `alter()`, `egoalterint()`, `tertius()` and `four` in the choice model and `ego()` in the rate model; `opportunities()` together with the attribute effects (`same`, `diff`, `sim`, `alter`, `tertius`, `egoalterint`) through the new option `oppframe(framename)`; `change()` in two-mode networks and together with `opportunities()`. What goldfish itself refuses stays refused, with a message.
- **Effects on other networks.** `netterms()`: `tie`, `trans`, `cycle`, `commonsender`, `commonreceiver` and `four` through another network, each with its own `window()`, `weighted` and `ignorerep`; the mixed closures `mixedtrans`, `mixedcycle`, `mixedcommonsender` and `mixedcommonreceiver` with the fitted network and an exogenous one, or with two exogenous ones; `tertiusdiff(var [net])`; `tie(netname)` accepts a network that changes over time; `ignorerep(effects)` for the built-in choice and choice_coordination effects.
- **Postestimation, as in goldfish.** `estat ic` (log likelihood, log likelihood per event, AIC, BIC), `estat eventlik` and `estat eventprob` (goldfish's `returnIntervalLogL` and `returnEventProbabilities`), `estat outliers` (Hampel, IQR, Top; goldfish's `examineOutliers`) and `estat changepoints` (PELT, AMOC and binary segmentation for mean and variance; goldfish's `examineChangepoints`). They run in the compiled plugin. The change point methods are a port of the R package changepoint 2.3 (GPL) and are credited in NOTICE. For the rate model with intercept the right-censored intervals that goldfish counts are reported (`llcens`, `r(ll_trailing)`, `r(n_cens)`) and AIC and BIC use goldfish's record count.
- After a fit its data stays in the frame `__nwdynam_fit`; the results carry `e(fitid)` and `e(estat_cmd)`, and `estat` stops with return code 498 when the results do not belong to the fit in memory.
- Estimates, standard errors and log likelihoods agree with goldfish 1.6.12 to 1e-6 on about 700 random and fixed event logs (coordination: 5e-4 and 10% for the standard errors, which are numerical in both programs); the postestimation values agree on 90 logs (per-event values to 1e-13, outlier and change point indices exactly). Plugin protocol 5: an older plugin is refused (r(3495)).

### Changed
- `opportunities()` is read as a string option; its error codes are unchanged.
- Fits that do not use the new features are not slower.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (full test suite and the extended tests). Windows and Linux: not yet tested for this version; the plugins are built by the build workflow, and until they are rebuilt those platforms stop with return code 3495.

## [1.7.0] - 2026-11-16

### Added
- Composition change, as in goldfish (the node attribute `present`): actors who join or leave the network over the log. Option `present(variable)` gives the presence at the start (0/1; default everybody), option `compchange(network)` the joins (value 1) and leaves (value 0) as an event network (one actor per event, sender = receiver). An absent actor is not in the risk set: it cannot send (rate sub-model, with and without `intercept`: it is left out of the hazard and its exposure), cannot be chosen (choice) and cannot be a member of a pair (choice_coordination); its statistics are still computed. A change takes effect for the events at its time and after it (a join at the time of an event lets the actor take part, a leave at that time excludes it); with `intercept` the waiting time that ends at an event is exposed to the actors present at its end, as goldfish. An event of an absent actor is an error (198), as in goldfish. Works for one- and two-mode networks, with `opportunities()` (intersection with the present actors; `present()` itself cannot be combined with `opportunities()`) and with `change()`.
- `e(ncompchange)`, `e(compchange)`, `e(present)` and a line in the output header when composition change is used.
- Estimates, standard errors and log likelihoods agree with goldfish 1.6.12 to 1e-6 on 67 random and fixed event logs (`dev/dynam_comp_crosscheck.R`, `dev/comp/`, test `cscripts/test_nwdynam_composition.do`).

### Changed
- The plugin protocol is 4: plugin binaries of 1.6.0 and earlier are recognised as too old (r(3495)). Until the build workflow has rebuilt the Linux and Windows plugins those platforms stop with r(3495). Fits without composition change give the same results as before to round-off and are not slower.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite; the goldfish cross-checks including 67 random logs with composition change). Windows and Linux: not yet tested; the plugins are built by the build workflow, and until they are rebuilt those platforms stop with return code 3495.

## [1.6.0] - 2026-11-16

### Added
- Time-varying covariates, as in goldfish (attribute change events linked to the nodes): option `change(variable network [variable network ...])` lets a covariate change over the log. The changes are an event network declared with `nwset` (one actor per event, sender = receiver, the new value as the event value); the dataset holds the values at the start. Option `increment` reads the event values as increments (goldfish: `increment` instead of `replace`). `same()`, `diff()`, `sim()`, `ego()`, `alter()`, `egoalterint()` and `tertius()` use the value that holds at the time of each event in the choice, rate (also with `intercept`) and choice_coordination sub-models. A change takes effect for the events after it; an event at exactly the time of a change still sees the old value (goldfish processes the events of a time point before the other event lists). With `intercept` the observation period runs from the first to the last event of either list, and every change splits the time between two events, as goldfish does.
- Estimates, standard errors and log likelihoods agree with goldfish 1.6.12 to 1e-6 on 42 random event logs with changes (the script and data are in `dev/`, the test is `cscripts/test_nwdynam_timevarying.do`).

### Changed
- The plugin protocol is 3: plugin binaries of 1.5.1 and earlier are recognised as too old (r(3495)). Until the build workflow has rebuilt the Linux and Windows plugins those platforms stop with r(3495). Fits without `change()` give the same results as before to round-off and are not slower.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite; the goldfish cross-checks including 42 random logs with changes). Windows and Linux: not yet tested; the plugins are built by the build workflow, and until they are rebuilt those platforms stop with return code 3495.

## [1.5.1] - 2026-11-16

### Faster
- Dense data (thousands of events among few actors): rate fits with `intercept` and choice fits with effects such as `trans` no longer spend up to a thousand likelihood evaluations after the maximum was already reached. `submodel(rate) intercept` on 5,000 events among 50 actors takes 0.02 s (was 2.4 s), `inertia trans` on 5,000 events among 100 actors 0.03 s (was 3.7 s). The estimates are the same to round-off.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite and, under Stata 17, the goldfish cross-checks). Windows and Linux: not yet tested; the plugins are those of 1.5.0.

## [1.5.0] - 2026-11-15

### Changed
- Every fit runs in the compiled plugin; the Mata evaluators are no longer part of the shipped fit path. Choice (all 17 effects with windows, weighted effects, two-mode networks, `opportunities()` and `tie()`), rate, rate with intercept (all effects, windows, weights, two-mode) now evaluate log likelihood, gradient and the exact Hessian natively (one generic kernel, one data block per fit, only theta travels per evaluation). The Hessian gives the standard errors; Newton-Raphson polishes the maximum of the concave models to round-off.
- A fit without `seed()` starts from a fixed point and is reproducible.
- When the plugin file is missing the command stops with return code 601; when the installed plugin is older than the package (or cannot be loaded) with return code 3495. Both name the file and say to reinstall. The warning line about the slower Mata engine is gone, with the engine.
- The plugin protocol is versioned (query code 99); plugin binaries of version 1.4.0 and earlier are recognised as too old. The Linux and Windows plugins are rebuilt by the build workflow; until then those platforms stop with r(3495).

### Faster
- Choice coordination on the manufacturing emails (undirected): `inertia` 5 s for 368 events (was 66 s) and 19 s for 1,019 events (was 253 s); `inertia trans` 9 s (was 113 s) and 31 s (was 449 s). The estimates are the same to the printed digits. Rate fits without `intercept` start from a Newton step: the default fit that took 59 s (about 2,000 evaluations) takes 0.8 s.

### Removed
- The Mata fit paths and evaluators. They live on as reference code of the test suite (`cscripts/ref_dynam_mata.do`) against which the plugin is compared on random cases.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite and, under Stata 17, the extended tests). Windows and Linux: not yet tested; the plugins are built by the build workflow.

## [1.4.0] - 2026-11-15

### Changed
- `nwdynam, submodel(rate) intercept` is fast. The first 14 days of the manufacturing emails (2,738 events, 129 actors) with `intercept indeg` take 0.2 seconds (before: no result after an hour); 50,000 events among 1,000 actors take a few seconds. All statistics of the real-time hazard are computed once instead of at every likelihood evaluation, and for `indeg` and `outdeg` (one-mode, no window, not weighted) the native plugin evaluates the likelihood without storing them.
- The estimate is found by Newton-Raphson on the analytic Hessian from a data-driven start (the log likelihood is concave), then polished until the gradient is at round-off level; one start instead of eight. Coefficients and standard errors are now the converged maximum-likelihood values to 12 printed digits. Fits made with earlier versions agreed with that maximum only to about 5 digits (the eight random starts ended at slightly different points), so a re-run can change the last digits of an old `intercept` result.

### Fixed
- Larger problems with `intercept` no longer fail: the optimizer used to leave the region where the hazard can be computed (`exp()` overflowed to missing and the hazard term was silently dropped), so no start converged and a Nelder-Mead search followed. The fit now stops a step at such a point and shortens it.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the nwdynam test suite, including the goldfish cross-checks under Stata 17) and on Linux (Ubuntu 22.04) with StataNow 19.5 MP, where the Linux plugin takes the native path and agrees with the Mata evaluator. Windows: not yet tested in Stata; the Windows plugin is built by the build workflow, and until it is tested there that platform may use the Mata evaluator, which is fast as well.

## [1.3.2] - 2026-11-13

### Changed
- Rebuilt for nwcommands 3.24.0; requires core 3.24.0 or later. No other change.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite), on Windows 11 with StataNow 19.5 MP and on Linux (Ubuntu 22.04) with StataNow 19.5 MP.

## [1.3.1] - 2026-11-13

### Changed
- Rebuilt for nwcommands 3.23.0; requires core 3.23.0 or later.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite) and on Linux (Ubuntu 22.04) with StataNow 19.5 MP. Windows: not yet tested in Stata.

## [1.3.0] - 2026-11-13

### Results that change (read before you re-run old analyses)
- Covariate effects on a dataset whose rows are not in actor order but carry `_nwnode`: the covariate is now matched to the right actors.
- Fits that silently used events with a missing time or a missing actor now stop with r(416).
- Valued event networks (a third variable in `nwset ..., eventtime()`): `weightedinertia`, `weightedrecip`, `weightedindeg` and `weightedoutdeg` sum the event values instead of counting events, as goldfish does. Values of 0 or below also change the unweighted effects: an event with the value 0 creates no tie, a negative value can remove one (every effect then sees the tie as absent), and the unweighted `indeg`/`outdeg` add the sign of each tie value (a tie below 0 counts -1), all as in goldfish. With only positive values the unweighted effects do not change; networks without a value variable do not change at all.
- Windows on a network with values of 0 or below use the summed value of the events within the window (a tie is present while it is above 0); with positive values nothing changes.

### Fixed
- Covariates (`same()`, `diff()`, `sim()`, `ego()`, `alter()`, `tertius()`, `egoalterint()`) are matched to the actors by their node labels (variable `_nwnode`, as after `nwload, xvars`) instead of by row position; a re-sorted dataset no longer misaligns them. Without `_nwnode` the data must still have one row per actor in actor order. An actor without a row, or with more than one, is refused (r(198)).
- A missing covariate value is refused (r(416)); before, `alter()`/`ego()` failed with r(111) after many optimizer lines and `same()` silently treated missing as a category of its own.
- Events with a missing time, or with a sender or receiver that is not an actor of the network (for example a missing id), are refused (r(416)); before, the native engines silently fitted them (an event without a time as the last event, a missing id as actor 0) and the other engines gave raw Mata errors.
- `opportunities()`: as in goldfish, the realised receiver of every event must be among the available actors (r(198)); before, the fit ran with a log likelihood of -1e300. Rows with a missing value are refused (r(416)) instead of a raw Mata error.
- A failed fit no longer leaves the user in the helper frame `__nwdynam_native` (which made every later `nwdynam` call fail with r(3598)); the helper scalar `__nwdynam_sequential` is removed after errors as well.
- A self-loop event gives its message and r(3498); the message was lost (Mata r(3300)) because it was too long for Mata's `_error()`.
- An effect with zero variance (for example a window that no earlier event falls into) is listed in `e(unidentified)` and in the note.
- The window options and `seed()` no longer treat -1 as "not given": windows must be positive numbers, `seed()` an integer from 0 to 2147483647 (r(198) otherwise).
- `nwdynam` without arguments (or with `level()` only) after an `nwdynam` fit redisplays the results, and `estimates replay` works (it used to refit the current network); new option `level()`, new stored result `e(header)`.
- Help: `seed()` describes which effect sets start from random values (the default effect sets use fixed starts).
- Event values are used as in goldfish (needs the core with event values, `nwset sender receiver value, eventtime(time)`): the value of an event is added to the tie (goldfish's `increment`); a tie exists while its value is above 0. Before, every event counted 1 and the values were ignored. Events with a missing value are refused (r(416)).
- A weighted effect can be combined with its own window (`weightedinertia inertiawindow()`, `weightedrecip recipwindow()`, `weightedindeg indegwindow()`, `weightedoutdeg outdegwindow()`): the summed value of the events within the window, as goldfish's `weighted = TRUE, window =`. Before, the combination was refused (r(198)).

### Changed
- Rebuilt for nwcommands 3.22.0; requires core 3.22.0 or later.

Mata library `lnwdynam.mlib` rebuilt with Stata 17 against the new core library (new functions `dynam_node_rows()` and `dynam_window_expire()`, self-loop message, event values). The native plugin `dynam_sim` takes the event values (used only when some value is 0 or below); rebuilt for macOS, Linux and Windows.

Tested on macOS with Stata 17 BE, Stata 18 BE and StataNow 19.5 MP (the full test suite). Windows and Linux: plugins built by CI; the Windows plugin was not yet tested in Stata for this release.

## [1.2.5] - 2026-11-12

### Fixed
- `opportunities()` combined with `same()`, `diff()`, `sim()`, `alter()`, `ego()`, `tertius()` or `egoalterint()` stops with a clear error (error 198). The combination was documented as not possible but was not checked.
- `tie()` and `opportunities()` on a two-mode network stop with a clear error (error 198), as the other effects not yet verified for two-mode networks do.
- `egoalterint()` with a string variable stops with a message that names the option and the variable (error 109).

## [1.2.4] - 2026-11-12

### Fixed
- **Help file checked against the code.** Without effect options the fit uses a fixed default set of effects (it does not fit every effect); the option table says which sub-models each effect applies to and which network types each sub-model takes; `seed()`, the native engine and the stored results are described as they are.
- **Help file layout.** Stored results wrap at the window width; very long source lines are broken up.
- **The help file describes the current behaviour only**: references to earlier versions and development notes are removed; the option `sequential` is described by what it does.

Results do not change.

## [1.2.3] - 2026-11-11

### Fixed
- **Long notes, warnings and messages wrap at word boundaries instead of being cut in the middle of a word.** Stata breaks a displayed line that is longer than the line size at exactly that column (for example in the error "nwdynam requires a network declared via nwset's eventtime() option ..." and the note on effects that are not identified); these messages are now SMCL paragraphs, which wrap between words. The wording is unchanged. Results do not change. Tested in `cscripts/test_nwdynam_msgwrap.do` at line size 80.

## [1.2.2] - 2026-11-09

Patch release. No change of the Mata library or the estimators.

### Changed
- The Windows plugin is compiled with optimisation (`/O2`), like those of the core, nwsaom and nwergm. The results are the same; the macOS and Linux plugins are unchanged. The full test suite passes on Windows.
- Affiliation in `nwdynam.pkg`: RWTH Aachen University.

### Added
- Documentation of `nwinstall, buildmata(nwdynam)` (nwcommands 3.14 and later): a Mata library that does not load in the user's Stata (built by a newer Stata, r(3499)) is rebuilt from the source in this repository with the Stata that runs. README and help file describe it. No change of the Mata library, the plugin or the estimators.

## [1.2.1] - 2026-10-31

Patch release. The build routine `nwinstall, buildplugin(nwdynam)` needs nwcommands 3.13.0; the package itself works with the core it worked with before. The compiled plugin is unchanged.

### Added
- The C source of the native plugin (`native/*.c`, `native/*.h`, the `Makefile`, StataCorp's plugin interface files and a build note `dynam_native_readme.txt`) is listed in `nwdynam.pkg` as ancillary files: `net get nwdynam, from(...)` copies it; `net install` still installs the compiled plugins only.
- If the shipped plugin does not work on a computer, `nwinstall, buildplugin(nwdynam)` (nwcommands 3.13.0 or later) fetches the source, compiles it there and installs the plugin (macOS and Linux).

## [1.2.0] - 2026-10-26

Minor release (documented behaviour changes). Works with nwcommands 3.6.0 and later; no change of the Mata library.

### Added
- Warning when the compiled plugin is missing: a fit of the original structural effect set (choice and rate sub-models) that falls back to the slower Mata engine now prints one warning line (once per call), as nwergm does. Results are unchanged. New helper `_nwdynam_slowwarn.ado`; new test `cscripts/test_nwdynam_fallbackwarn.do`.

### Changed
- Help text: the Mata fallback is no longer described as transparent; it states the warning and which fits can warn.

## [1.1.4] - 2026-10-26

Patch release. Works with nwcommands 3.6.0 and later; no change of the Mata library.

### Changed
- Help text corrected: the long list of internal verification scripts is replaced by a short statement of the comparison with goldfish; the README install lines start with nwinstall. No change in code or results; the Mata library is unchanged.

## [1.1.3] - 2026-10-26

Patch release. Works with nwcommands 3.6.0 and later; no change of the Mata library.

### Results that change (read before you re-run old analyses)
- `tie()` counts a tie only when its value is greater than 0, as goldfish's I(x_ij>0); negative values used to count as ties. Covariate networks without negative values give the same results as before.

### Added
- `cscripts/test_coreaudit.do`: valued event weights and `tie()` under core 3.6.0.

## [1.1.2] - 2026-10-26

### Changed
- Rebuilt for nwcommands (core) 3.6.0: the Mata library is compiled against the core library of 3.6.0 (built with Stata 17) and the package needs nwcommands 3.6.0 or later. No change in results.
- Install lines and the README: `nwinstall, everything` installs the core and the four model packages in one step. Supported Stata versions: 17 or later, tested with Stata 17 and 18.

## [1.1.1] - 2026-10-25

### Fixed
- The Mata library was built with Stata 18 and could not be used in Stata 17; it is now built with Stata 17 (works in Stata 17 and 18). Needs nwcommands 3.5.1 or later.

### Changed
- Requirements stated plainly: Stata 17 or later; tested with Stata 17 and Stata 18.


## [1.1.0] - 2026-10-25

### Results that change (read before you re-run old fits)
- Events with the same time are now simultaneous, as in goldfish: the statistics of all events at a time point come from the state before that time point. Estimates for data with tied times change; the option `sequential` keeps the old way. Data without tied times are unchanged.

### Changed (results change for data with tied event times!)
- Events with an identical time are now simultaneous, as in goldfish 1.6.12: the statistics of all events of a time point are computed from the state before that time point (all effects, windows and the rate intercept included, and the coordination model), and the state is updated after the whole group. Before, tied events were processed one after another in input order. Example: `socialevolution` with `ceil(t/20)`: log likelihood -652.83 and inertia 5.00 before, -1595.33 and 3.54 now (goldfish: -1595.33, 3.54). Data without ties give identical results (byte-identical estimates on the book examples).
- New option `sequential` keeps the old behaviour. New `e(ntimes)`, `e(ties)`; a note under the header says how many events share a time point.
- Native plugin: new codes 11/12 (simultaneous events, extra variable v4). Linux and Windows plugins must be rebuilt by CI before release; an older plugin stops with an error on tied data instead of giving wrong numbers.

### Fixed
- Events with the same time keep the order in which the network lists them. Before, the events were sorted with `order()` on the time column alone, which breaks ties by the state of the random-number generator, so the same data with tied times gave a different tie history and a different likelihood from run to run.

## [1.0.4] - 2026-10-25

### Changed
- Rebuilt for nwcommands 3.4.0; requires core 3.4.0 or later.
- A fit with effects that are not identified (no standard error, or a standard error above 100, as in a very short event log) now prints a note under the table and returns the effect names in `e(unidentified)`.

### Fixed
- `nwdynam version` as the first command of a fresh session no longer reports the native plugin as "missing": the Mata library is indexed before the plugin path is looked up.

## [1.0.3] - 2026-10-25

### Fixed
- Works right after a fresh install: the Mata library is indexed automatically. Before, the first call in a new Stata session could stop with r(3499) "... not found" (or, for nwergm, "class ErgmGraph undefined") until `mata mlib index` was run by hand.

## [1.0.2] - 2026-10-25

### Changed
- The version is shown at the bottom of each help file.

## [1.0.1] - 2026-10-24

### Changed
- The help file documents the `nwdynam version` subcommand and the stored results `r(version)`, `r(date)`, `r(core)`, `r(license)` and `e(version)`.

## [1.0.0] - 2026-10-24

First release of nwdynam as a package of its own.

### Added
- Dynamic network actor models (DyNAM; Stadtfeld and Block 2017) for event-time networks declared with `nwset ..., eventtime()`: rate and choice submodels.
- Network effects (inertia, reciprocity, in- and out-degree, transitivity, cycles, common sender/receiver, four-cycles, node-level transitivity), covariate effects (same, diff, sim, ego, alter, tertius, ego-alter interaction), windowed and weighted effects, opportunity sets, ties as effects, intercepts; two-mode networks.
- Native C plugin (macOS, Linux, Windows) for the likelihood and gradient.
- `nwdynam version`: package, core and plugin versions; `e(version)` in the results.
- Compared with and validated against R goldfish 1.6.12.
