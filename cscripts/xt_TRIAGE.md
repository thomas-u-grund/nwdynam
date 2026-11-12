# nwdynam extended tests: triage (2026-11-12, Stata 17 BE, nwdynam 1.2.5 + core of the same date)

Repro do-files and logs: `$S/xtrun/models/repro/nwdynam/` (`dyn_repro.do` 8-event toy, `dyn_repro2.do` W40,
`dyn_repro3.do` frame leak and same(), `tm_str.do`/`tm_str2.do` two-mode string ids, `gf.R` goldfish recheck).
Test logs: `$S/xtrun/models/nwdynam_xt_dyn_*.log`. Owner decisions of 2026-11-12 are recorded below.

## Suspects

| Id | Verdict | Severity | One-line repro | Cause |
|---|---|---|---|---|
| S3 | CONFIRMED | wrong results silently | W40 plus a row with missing sender (or missing receiver): `nwdynam w, seed(1)` and `submodel(rate)` fit 41 events, no message; the event list holds actor 0 | core build_eventlist (unw_core.do:484) gives index 0; the native unit1 engines do not check. Multi engine (`inertia trans`): raw r(3301). On 8 events, missing receiver: r(111) after "flat region" lines (and then I11) |
| S4 | CONFIRMED | wrong results silently | `replace t = . in 10`: choice fit uses the timeless event as the last one (e(N) = 40); `submodel(rate) intercept`: r(111) after 8 x "initial values not feasible" | core nwset.ado:322-323 accepts missing times; DynamState::init sorts them last (unw_dynam.do:180) |
| S5 | CONFIRMED | wrong results (visible: ll = -1.0e+300) | opportunities() without the realised receiver of event 1: fit runs, ll -1e300, inertia -13 (toy) | nwdynam.ado:692-701 checks ranges only; goldfish refuses ("Active node 3 not available in event 1") |
| S12 | CONFIRMED | wrong results silently (goldfish deviation) | valued event network `nwset s r w, eventtime(t)`, `nwdynam wval, weightedinertia`: 0.2889 / ll -37.7174 (= unit increments); goldfish with increment = w: 0.1389639466 / -39.2302174436 (rechecked in R today) | weighted effects count events and ignore the values; the help says "matching goldfish weighted(TRUE)" and does not document the deviation |
| S13 | CONFIRMED (owner 2026-11-12: match covariates via `_nwnode` when present) | wrong results silently | `nwload, xvars`, `gsort -_nwnode`, `nwdynam w40, alter(x4)`: -0.357 instead of -0.496 | nwdynam.ado:603-640 read by row position |
| I11 | CONFIRMED (notes item 11, reproduced without fault injection) | session broken, data appear lost | 8 events + one with a missing receiver: `nwdynam c2` -> r(111); the user is left in frame `__nwdynam_native` (4 helper variables, own data not visible); every later `nwdynam` call -> "frame __nwdynam_native already defined", r(3598), until `frame change default; frame drop __nwdynam_native` | unw_dynam.do:457-460 create the frame, 519-520 / 683-684 restore it only on the normal path |
| S1 | CONFIRMED | error only (raw Mata error) | ids >= 10,000,000 (`10000000 + id`): `_error(): 3300 argument out of range` | core nwset.ado:285-286 `strofreal()` %9.0g vs %18.0g node names (nwfromedge.ado:197, 246); index 0/0 hits the self-loop guard -> N2 |
| S2 | CONFIRMED | error only | ids 1e9..4e9: same r(3300) | as S1 |
| S6 | CONFIRMED | error only (raw Mata error) | opportunities() data with a missing value: `<istmt>: 3200 conformability error` | nwdynam.ado:694, 699 `if !missing()` skip them; loop at 704 subscripts with missing |
| S11 | CONFIRMED (owner: covariates must be checked) | alter()/ego(): error only (misleading); same(): wrong results silently | x with one missing: `alter(x)` -> r(111) after "initial values not feasible"/"simplex delta required"; `same(x)` runs and treats missing as its own value (same b as with 9 in its place) | nwdynam.ado:603-640 no missing check |
| N2 (new) | CONFIRMED | error only (message lost) | self-loop event: `nwdynam loop` (all sub-models) -> `_error(): 3300 argument out of range`, no message | unw_dynam.do:167: Mata `_error(errtxt)` fails with r(3300) on texts longer than about 110 characters. Same bug in nwrem |
| N3 (new, core) | CONFIRMED | error only here; core network inconsistent | `dev/dynam_crosscheck_twomode_events.csv` (senders "P..", receivers "O..") `nwset sender receiver, twomode eventtime(time)` -> `nwdynam` r(198) "60 event(s) violating" mode direction | core nwset twomode with string ids sorts node names across both modes (O1 O2 P1 P2 P3) but keeps Level 1 = number of senders; static `nwset ..., twomode` shows the same mismatch (nwtomata missing blocks follow 2+3, metadata says 3+2). Numeric codes, or labels where senders sort first, work. Check nwergm bipartite terms on such networks |
| S7 | CONFIRMED | cosmetic | `nwdynam one` (1 event) -> r(2001); scalar `__nwdynam_sequential` remains | nwdynam.ado:657 created before the checks, dropped only at 816 |
| S8 | CONFIRMED | cosmetic (reporting) | `nwdynam w40, inertiawindow(0.5)`: V = 0, se shown ".", e(unidentified) empty | nwdynam.ado:855-860 test `missing(se) \| se > 100` only (nwrem also tests v <= 0) |
| S9 | CONFIRMED | minor | `inertiawindow(-1)` (also recip/indeg/outdegwindow(-1)) = no window, silently; -2 refused | nwdynam.ado:16, 362-422 sentinel -1 |
| S10 | CONFIRMED | minor | `set seed 5` / `set seed 6`, `nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(-1)`: e(b) differ | nwdynam.ado:16, 642-644 sentinel -1 |
| 12 (help) | CONFIRMED | cosmetic (help text) | default fit with `set seed 1` and `set seed 999`: identical e(b) | nwdynam.sthlp:57 says most effect sets start from random values; the default (unit1) engines use fixed starts |
| 13 | NOT A BUG (dead check reachable in practice) | - | the check at nwdynam.ado:765-772 is reachable through N3 (string ids); with `nwset o p, twomode` the roles just swap | keep the check |
| 14 (core) | CONFIRMED | cosmetic (message) | `nwdynam loop g` -> "Network loop g not found", r(482) | core _nwsyntax.ado:70 |
| N1 (new) | DECIDED: fix (owner 2026-11-12: add estimates replay to all four packages, low priority) | error only / surprising | `estimates replay` runs `nwdynam` without arguments, which REFITS the current network with the default effects (or fails, e.g. r(198) when the current network is undirected) | no replay branch in nwdynam.ado |

Counts: CONFIRMED 19 (S1-S13, I11, N2, N3, item 12, item 14; S13 by owner decision), NOT A BUG 1 (item 13), DECIDED: fix 1 (N1, replay).

## Test-code fixes

- xt_dyn_edges.do: `nwdynam w40, sequential` after w40 had been cleared -> redeclare w40; id map `10000000 + 2*@ - 1`
  recorded as CONFIRMED BUG S1; the two-mode reference fit now uses numeric codes, the string P../O.. fit is recorded
  as CONFIRMED BUG N3.
- xt_dyn_estimation.do: `e(effects)` includes `Intercept` (help: e(b) columns named per e(effects)) -> expectation
  fixed; `estimates replay` non-fatal (N1); matrix literal in `mreldif(xt_keep_matrix, (1, 2))` -> stored matrix;
  two-mode CSV via numeric codes (N3).
- xt_dyn_errors.do: second `nwset` after the first had cleared the data (twice) -> reload; two-mode CSV via numeric
  codes (N3); the five r(3498) self-loop checks recorded as CONFIRMED BUG N2.
- xt_dyn_suspected.do: verdict markers; S11 now asks for a clean refusal (198/416) for alter, ego and same; frame
  cleanup after every check; new checks N2 and I11.

## Remaining failures per file

| File | Result |
|---|---|
| xt_dyn_edges.do | ALL PASS apart from CONFIRMED BUG S1, N3 |
| xt_dyn_closedform.do | ALL PASS |
| xt_dyn_clogit.do | ALL PASS |
| xt_dyn_estimation.do | ALL PASS (N1 replay noted; owner decided to add it) |
| xt_dyn_options.do | ALL PASS |
| xt_dyn_errors.do | ALL PASS apart from CONFIRMED BUG N2 |
| xt_dyn_suspected.do | ends r(9) by design: S1-S13, N2, I11 all fail (all confirmed) |

## Fixes (2026-11-12, Unreleased)

Fixed in nwdynam: S3, S4 (r(416)), S5 (r(198), as goldfish), S6 (r(416)), S7, S8, S9, S10, S11 (r(416)), S13 (match by
`_nwnode`), I11 (frame restored after any failed fit), N2, item 12 (help), N1 (replay, level()). S12: documented as a
deviation plus a note after the fit; using the values as goldfish does needs per-event values that the core event list
does not store (core change). Core-caused, left for the core: S1, S2 (nwdynam refuses the unmatched events cleanly,
r(416)), N3 (two-mode string ids mode split), item 14 (_nwsyntax message). All xt_dyn_*.do end with ALL PASS under
Stata 17.
