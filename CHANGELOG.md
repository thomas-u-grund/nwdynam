# Changelog

All notable changes to nwdynam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- `nwdynam version` as the first command of a fresh session no longer reports the native plugin as "missing": the Mata library is indexed before the plugin path is looked up.

### Changed
- A fit with effects that are not identified (no standard error, or a standard error above 100, as in a very short event log) now prints a note under the table and returns the effect names in `e(unidentified)`.

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
