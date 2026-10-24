# nwdynam

`nwdynam` fits dynamic network actor models (DyNAM) to relational event data declared as an event-type temporal network with `nwset ..., eventtime()`: a rate model (when does an actor act) and a choice model (whom does the actor choose), with network, covariate, window and coordination effects. The hot likelihood/gradient evaluation runs in a native plugin.

Part of the nwcommands family of Stata packages for network analysis. nwdynam is a separate package with its
own repository; it needs the core package nwcommands but no other model package.

## Install

1. Install the core package first (Stata): `net install nwcommands, from(<URL of the nwcommands distribution>)`
   then run `nwinstall`. Required: nwcommands (core), distribution of 5sep2026 or later.
2. Install this package:

       net install nwdynam, from(https://raw.githubusercontent.com/<account>/nwdynam/main/)

   (`<account>` is a placeholder for the GitHub account that hosts the repository.)

## Requirements

* Stata 14 or later (StataBE is enough).
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

## Native plugin

The package ships a native Stata plugin (dynam_sim.c) for macOS (universal), Linux and Windows (built by the GitHub Actions workflow `.github/workflows/build-plugins.yml`; the macOS build can be reproduced with `make -C native macos`). Platform without a plugin: see the help file.

## Rebuilding the Mata library

`lib/build.do` rebuilds `lib/lnwdynam.mlib` from `unw_dynam.do` (run it from the repository root: `StataBE -e do lib/build.do`).

## License

GNU General Public License version 3 (see `LICENSE`); see `NOTICE.md` for what the package is derived from and
how to cite it.
