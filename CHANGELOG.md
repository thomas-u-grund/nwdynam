# Changelog

All notable changes to nwdynam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

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
