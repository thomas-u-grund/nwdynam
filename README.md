# nwdynam

<!--version-->**Version 1.11.2** (released 26 November 2026)<!--/version-->

`nwdynam` fits dynamic network actor models (DyNAM) to relational event data declared as an event-type temporal network with `nwset ..., eventtime()`: a rate model (when does an actor act) and a choice model (whom does the actor choose), with network, covariate, window and coordination effects. The hot likelihood/gradient evaluation runs in a native plugin.

Part of the nwcommands family of Stata packages for network analysis. nwdynam is a separate package with its
own repository; it needs the core package nwcommands but no other model package.

## Install

1. Install the core package first (Stata): `net install nwcommands, from(https://nwcommands.org/net/) replace`
   (alternative: `from(https://raw.githubusercontent.com/thomas-u-grund/nwcommands/main/)`). Required: nwcommands (core), version 3.29.0 or later.
2. Install this package from the core: `nwinstall, modelpkg(nwdynam)`.
   `nwinstall, everything` installs the core and all four model packages in one step; `nwinstall, all` installs the core only and `nwinstall, models` adds the model packages to it.

`nwinstall` installs the latest release (the state of the main branch). Without `nwinstall`, the package can be installed directly from GitHub:

       net install nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main/) replace

To install a specific version, point `from()` at the release tag (GitHub serves every tag under its own raw URL), for example version 1.1.3:

       net install nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/v1.11.2/) replace

`net describe nwdynam` shows the version and distribution date, `nwdynam version` the package, core and plugin versions, and `e(version)`
holds the package version in the stored results. Release history: `CHANGELOG.md`; how releases are made: `RELEASING.md`.


## Requirements

* Stata 17 or later. Tested on macOS with Stata 17, 18 and 19 (Stata/BE is enough); the Windows and Linux plugins are built in the same way. The Mata libraries are built with Stata 17, the oldest tested version, because a library built by a newer Stata is ignored by older ones; Stata 16 and older have not been tested and are not expected to work.
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

## Time-varying covariates

Node covariates may change while the events happen, as in goldfish (attribute change events linked with `linkEvents()`). The dataset holds
the values at the start of the log; the changes are an event network of their own (one actor per event, the new value as the event value),
passed with `change(variable network)`:

    nwset id id2 newvalue, eventtime(time) selfloop name(xchange)     // id2 = copy of id
    nwdynam net, inertia alter(x) change(x xchange)

The effects `same`, `diff`, `sim`, `ego`, `alter`, `egoalterint` and `tertius` then use the value that holds at the time of each event, in the
choice, rate (also with `intercept`) and choice_coordination sub-models; `increment` reads the values of the change events as increments. See
`help nwdynam`. The R script and data of the comparison with goldfish are `dev/dynam_tvc_crosscheck.R` and `dev/tvc/`.

## Composition change

Actors may join or leave the network while the events happen, as in goldfish (the node attribute `present`). An absent actor is not in the risk
set (it cannot send in the rate sub-model, cannot be chosen in the choice sub-model, cannot be in a pair of choice_coordination); its statistics
are still computed. `present(variable)` gives the presence at the start (0/1, one row per actor; default: everybody is present),
`compchange(network)` the joins (value 1) and leaves (value 0) as an event network with one actor per event:

    nwset id id2 joined, eventtime(time) selfloop name(membership)     // id2 = copy of id
    nwdynam net, inertia indeg present(member0) compchange(membership)

A change counts for the events at its time and after it; an event of an absent actor is an error, as in goldfish. Estimates, standard errors and log
likelihoods agree with goldfish 1.6.12 to 1e-6 on 67 random and fixed event logs (choice, rate, rate with intercept, choice_coordination, one- and
two-mode, together with time-varying covariates). See `help nwdynam`. The R script and data are `dev/dynam_comp_crosscheck.R` and `dev/comp/`.

Windowed rate effects (`indegwindow()`, `outdegwindow()`) work together with `intercept`; `tertiuswindow()` puts a window on `tertius()`; `choice_coordination` accepts `inertiawindow()`, `indegwindow()`, `weightedinertia` and `weightedindeg`; two-mode networks accept `alter()`, `ego()`, `egoalterint()`, `tertius()` and `four`; `opportunities()` can be used with covariate effects when its variables are in another frame (`oppframe()`); `change()` works in two-mode networks and with `opportunities()`. These agree with goldfish 1.6.12 on 306 random event logs (`dev/dynam_win_crosscheck.R`, data in `dev/win/`); see `help nwdynam` for what is not available.

## Effects on other networks

`netterms()` adds the effects of goldfish that read other networks, with their own window: `tie`, `trans`, `cycle`, `commonsender`, `commonreceiver`, `four`, the mixed closures `mixedtrans`, `mixedcycle`, `mixedcommonsender`, `mixedcommonreceiver` (two networks, each the fitted network, an ordinary network or a network that changes over time with its events), and `tertiusdiff`. Options per term: `window(#)`, `weighted` (`tie`), `ignorerep`; `ignorerep(effects)` does the same for the built-in effects. `tie(netname)` also accepts a network that changes over time.

    nwdynam net, inertia netterms(tie(advice, weighted window(30)) mixedtrans(net advice) tertiusdiff(age))
    nwdynam net, inertia trans ignorerep(trans)

Estimates, standard errors and log likelihoods agree with goldfish 1.6.12 to 1e-6 on about 300 random event logs (`dev/dynam_ex_crosscheck.R`, data in `dev/ex/`, test `cscripts/test_nwdynam_exogenous.do`); see `help nwdynam` for what is not available.

## Postestimation

After `nwdynam`: `estat ic` (log likelihood, AIC, BIC), `estat eventlik` and `estat eventprob` (log likelihood and probabilities of every event, goldfish's `returnIntervalLogL` and `returnEventProbabilities`), `estat outliers` (Hampel, IQR, Top rules of `examineOutliers`) and `estat changepoints` (PELT, AMOC, binary segmentation for the mean or the variance, as `examineChangepoints`; the algorithms are a port of the R package changepoint, see `NOTICE.md`). Everything runs in the native plugin. Per-event values and probabilities agree with goldfish to 1e-13 on 90 random event logs (`dev/dynam_post_crosscheck.R`, test `cscripts/test_nwdynam_postestimation.do`); see `help nwdynam_estat`.

## Native plugin

The package ships a native Stata plugin (dynam_sim.c) for macOS (universal), Linux and Windows (built by the GitHub Actions workflow `.github/workflows/build-plugins.yml`; the macOS build can be reproduced with `make -C native macos`). The Windows plugin passes the full test suite (Stata 19 MP, 9 November 2026); the Linux plugin has not been tested on a Linux machine. Every fit needs the plugin; a platform without one stops with r(601) (see the help file). The C source is part of the package (`net get nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main)` copies it into the current directory, with `dynam_native_readme.txt` explaining the build) and is in `native/` of this repository. If the shipped plugin does not work on your computer (missing for the platform, older than the `.ado` files, does not load), `nwinstall, buildplugin(nwdynam)` compiles it there (macOS and Linux; needs a C compiler and `make`; restart Stata afterwards). Plugin code 11/12 (simultaneous events) needs a plugin built from this version; the Linux and Windows binaries are rebuilt by the workflow.

## Rebuilding the Mata library

`lib/build.do` rebuilds `lib/lnwdynam.mlib` from `unw_dynam.do` (run it from the repository root: `StataBE -e do lib/build.do`).

## Manual

A user manual (PDF in `manual/`, with the LaTeX source and the do-files of its examples) is archived at Zenodo: <https://doi.org/10.5281/zenodo.24398099> (all versions). Cite the manual by that address together with the package version printed by `nwdynam version`.

## License

GNU General Public License version 3 (see `LICENSE`); see `NOTICE.md` for what the package is derived from and
how to cite it.

## If the Mata library does not load

The Mata library (lnwdynam.mlib) is built with Stata 17. A library built by an older Stata loads in every later Stata; one built by a newer Stata is ignored by an older one (r(3499), "too new to be run by this version of Stata"). Build it with your own Stata: `nwinstall, buildmata(nwdynam)` (or `nwinstall, buildmata(all)` for the core and every installed package). The command fetches the Mata source from this repository, compiles it, installs the library and keeps the old file as `.prev`. It clears the data and the Mata objects in memory, so save first; restart Stata afterwards. Needs nwcommands 3.14 or later.
