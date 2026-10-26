# NOTICE

<!--version-->**Version 1.1.3** (released 26 October 2026)<!--/version-->

nwdynam is Copyright (C) 2026 Thomas U. Grund and is distributed under the GNU General Public License, version 3 (file `LICENSE`).

## What this package is derived from

`nwdynam` implements the DyNAM of Stadtfeld and Block (2017) and was developed from the R package `goldfish` (effect definitions, likelihoods and estimation), whose license is GPL (>= 3); `nwdynam` is therefore distributed under the GNU General Public License version 3.
* **goldfish** (R package, version compared: 1.6.12). Authors: James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors. License: GPL (>= 3).
  https://cran.r-project.org/package=goldfish
  Reference: Stadtfeld, C. and Block, P. (2017). Interactions, actors, and time: Dynamic network actor models for relational events. Sociological Science 4, 318-352. doi:10.15195/v4.a14

## Validation

Coefficients and log-likelihoods were compared with `goldfish` (version 1.6.12) on the bundled Social_Evolution data and on small hand-traceable examples; the cross-check scripts are in `dev/` (`dev/dynam_unit*_crosscheck.R` and the matching `.do` files) and the regression tests in `cscripts/test_nwdynam_*.do`. The native plugin is certified against the Mata reference implementation to floating-point round-off.

## AI assistance

The software was developed with Claude Code (Anthropic), an AI coding assistant, working under the author's direction. It helped port algorithms from goldfish into Stata, Mata and C and helped write the tests. The author specified the scope, reviewed the code, and is responsible for the package. Its correctness rests on validation, not on the assistant: the results are compared numerically with goldfish, and the test suite reruns these comparisons.

## How to cite

Grund, T. U. (2026). nwdynam: a Stata package (version 1.1.3). GitHub. https://github.com/thomas-u-grund/nwdynam
Please also cite the core package (Grund, T. U., nwcommands) and the work this package is derived from (see above).
