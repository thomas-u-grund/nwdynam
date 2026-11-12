# Changelog

All notable changes to nwdynam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

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
- Help: `seed()` describes which effect sets start from random values (the default effect sets use fixed starts); weighted effects on a valued event network count events, a documented deviation from goldfish (which adds the event values); a note says so after the fit.

### Results that change
- Covariate effects on a dataset whose rows are not in actor order but carry `_nwnode`: the covariate is now matched to the right actors.
- Fits that silently used events with a missing time or a missing actor now stop with r(416).

Mata library `lnwdynam.mlib` rebuilt with Stata 17 (new function `dynam_node_rows()`, self-loop message). No plugin change.

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
