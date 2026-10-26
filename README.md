# nwdynam

<!--version-->**Version 1.1.3** (released 26 October 2026)<!--/version-->

`nwdynam` fits dynamic network actor models (DyNAM) to relational event data declared as an event-type temporal network with `nwset ..., eventtime()`: a rate model (when does an actor act) and a choice model (whom does the actor choose), with network, covariate, window and coordination effects. The hot likelihood/gradient evaluation runs in a native plugin.

Part of the nwcommands family of Stata packages for network analysis. nwdynam is a separate package with its
own repository; it needs the core package nwcommands but no other model package.

## Install

1. Install the core package first (Stata): `net install nwcommands, from(<URL of the nwcommands distribution>)`
   then run `nwinstall, everything` (the core and all four model packages in one step) or `nwinstall, all` (core only). Required: nwcommands (core), version 3.6.0 or later.
2. Install this package:

       net install nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main/) replace

`net install nwdynam` installs the latest release (the state of the main branch). To install a specific version, point `from()` at the
release tag (GitHub serves every tag under its own raw URL), for example version 1.0.0:

       net install nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/v1.1.3/) replace

`net describe nwdynam` shows the version and distribution date, `nwdynam version` the package, core and plugin versions, and `e(version)`
holds the package version in the stored results. Release history: `CHANGELOG.md`; how releases are made: `RELEASING.md`.


## Requirements

* Stata 17 or later. Tested with Stata 17 and Stata 18 (Stata/BE is enough). The Mata libraries are built with Stata 17, the oldest tested version, because a library built by a newer Stata is ignored by older ones; Stata 16 and older have not been tested and are not expected to work.
* The core nwcommands package (see above), including its Mata library `lnwcommands.mlib`.
* Platform plugins are shipped for macOS, Linux and Windows.

## Tests

The tests run against an installed or checked-out core. Point the environment variable NWCORE to the
directory of the nwcommands (core) package (the directory that contains `unw_core.do` and the core `.ado`
files) and run, from the root of this repository:

    NWCORE=/path/to/nwcommands cscripts/run_tests.sh

Each test is a do-file in `cscripts/` that prints `PASS` or `ALL OK` lines; the runner prints the number of PASS
lines and any Stata error per test. Stata is run with `-e`. Validation scripts that compare with the R
reference implementation are in `dev/`.

## Events with identical times

Events with exactly the same time are simultaneous, as in goldfish: the statistics of all events of a time point come from the state before
that time point, and the state is updated after the whole group. The order of tied rows does not matter. Option `sequential` restores the
behaviour before 1.1.0 (tied events one after another). Data without ties give the same results as before.

## Native plugin

The package ships a native Stata plugin (dynam_sim.c) for macOS (universal), Linux and Windows (built by the GitHub Actions workflow `.github/workflows/build-plugins.yml`; the macOS build can be reproduced with `make -C native macos`). Platform without a plugin: see the help file. Plugin code 11/12 (simultaneous events) needs a plugin built from this version; the Linux and Windows binaries are rebuilt by the workflow.

## Rebuilding the Mata library

`lib/build.do` rebuilds `lib/lnwdynam.mlib` from `unw_dynam.do` (run it from the repository root: `StataBE -e do lib/build.do`).

## License

GNU General Public License version 3 (see `LICENSE`); see `NOTICE.md` for what the package is derived from and
how to cite it.
