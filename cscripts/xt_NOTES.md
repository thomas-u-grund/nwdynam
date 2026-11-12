# Extended nwdynam tests (xt_*.do): notes

DRAFTS, NOT YET RUN. Written 2026-11-12 against nwdynam 1.2.5 (main branch of the nwdynam repo, worktree branch
tests-extended) and the core at the same date, without starting Stata. A first run will probably show typing errors in
the test code itself. Fix those before you read any failure as a finding.

Run the files from the repository root with the usual runner, which sets the adopath to the core, this repository and
lib/plugins/<os>. The files contain no sysdir/adopath lines. They load `unw_dynam.do`, so they test the source and
not the mlib, and they use `dev/dynam_crosscheck_twomode_events.csv`. Each file ends with `di "ALL PASS"`, except
`xt_dyn_suspected.do`: if any suspected problem is confirmed, it lists the failing checks and stops with `exit 9`.

## Files

| File | Runtime (Stata 17 BE, estimate) | Coverage |
|---|---|---|
| xt_dyn_edges.do | ~40 s | Unsorted rows; time 0, negative times and monotone transforms (no window, no intercept); intercept with shifted times; tied times ceil(t/4) for choice, rate and intercept (goldfish); rows shuffled within a time point; sequential = distinct times; all events at one time point; string, offset and permuted ids; a valued event network; two-mode label variants; undirected coordination with pairs listed as (a,b) or (b,a), with and without ties. |
| xt_dyn_closedform.do | ~20 s | Hand scores and lls of tiny logs, sequential and simultaneous (choice, rate, coordination; multi engines = unit1 engines); binomial closed forms: rate ego(x01), intercept with ego, choice alter(x4), choice tie(nxt), the same with tied times; smallest data (2 events, 2 actors) and 60 actors. |
| xt_dyn_clogit.do | ~90 s | Own Mata builder for the choice and rate statistics (structural, closure, weighted, windowed, same/diff/alter/egoalterint/tertius/sim, nodetrans) with clogit on n = 6, N = 200 (rseed 20261113); distinct and tied times (simultaneous updates); sequential = distinct; opportunities() as dropped candidate rows; actor-label permutation. |
| xt_dyn_estimation.do | ~60 s | Seed: native engines independent, multi engine reproducible, seeds 1-5 reach the same optimum (= goldfish, 8 effects); every documented e() entry, including titles, e(ties) none/simultaneous/sequential and e(ntimes); `nwdynam version` r(); postestimation on choice native, choice multi, rate, rate+intercept, coordination, two-mode choice and rate; after a fit the data, c(frame), the event list, user scalars and matrices are unchanged, with no `__nwdynam_native` frame and no leftover scalars. |
| xt_dyn_options.do | ~60 s | Windows: inclusive boundary (5, 5.5, 4.5, 4.999999 against goldfish), window >= log span = no window, times x3 with window x3, time shift; weighted effects (goldfish); opportunities(): full sets = none, without sender rows, risk set of one (ll + ln 3), random row order, other names and types; intercept with times x2 (Intercept - ln 2, ll - 40 ln 2); option order; 8-effect goldfish reference; all 17 choice effects, all 5 rate effects, intercept column first. |
| xt_dyn_errors.do | ~20 s | Exact r() for every error exit of nwdynam.ado (198 list in the header, 109, 102, 103, 111, 482, 2001, 3498). After each error, e(b), e(ll), e(cmd) and the data are unchanged. |
| xt_dyn_suspected.do | ~30 s | S1-S13, the suspected bugs below. |

Total: about 5.5 minutes.

The clogit builder in xt_dyn_clogit.do follows the help-file definitions. Two of them are vague in the help file and
were taken from the implementation: `nodetrans` (number of transitive triads in which the actor is the source:
rowsum(A :* A*A)) and `tertius` (isolates get 0, as goldfish 1.6.12 actually does; see unw_dynam.do:969-989). If one
of those blocks fails, check the definition before you suspect the engine.

## W40 reference sequence

These are 40 events at t = 1..40 on 4 actors. They are generated in R by `set.seed(42); s <- sample(1:4, 40, TRUE);
r <- sapply(s, function(x) if (runif(1) < .5) setdiff(1:4,x)[1 + (x %% 3)] else sample(setdiff(1:4,x), 1))`:
`1 3;1 3;1 3;1 3;2 1;4 2;2 4;2 4;1 3;4 3;3 1;4 2;3 2;4 2;1 3;1 3;2 4;4 2;2 4;2 4;3 1;3 1;1 2;1 3;3 4;4 1;1 3;3 1;1 2;1 3;4 2;2 1;4 2;4 2;3 2;2 3;1 3;2 4;4 2;2 3`

## R code for the references (goldfish 1.6.12)

```r
suppressMessages(library(goldfish))
ev0 <- read.table(text=gsub(";", "\n", "1 3 1;1 3 2;1 3 3;1 3 4;2 1 5;4 2 6;2 4 7;2 4 8;1 3 9;4 3 10;3 1 11;4 2 12;3 2 13;4 2 14;1 3 15;1 3 16;2 4 17;4 2 18;2 4 19;2 4 20;3 1 21;3 1 22;1 2 23;1 3 24;3 4 25;4 1 26;1 3 27;3 1 28;1 2 29;1 3 30;4 2 31;2 1 32;4 2 33;4 2 34;3 2 35;2 3 36;1 3 37;2 4 38;4 2 39;2 3 40"))
n <- 4; lab <- paste0("A",1:n)
nodes <- defineNodes(data.frame(label=lab, present=TRUE, stringsAsFactors=FALSE))
fit <- function(tt, form, sub, directed=TRUE, inc=rep(1,40), init=list()) {
  ev <- data.frame(time=tt, sender=lab[ev0[,1]], receiver=lab[ev0[,2]], increment=inc, stringsAsFactors=FALSE)
  net <- defineNetwork(nodes=nodes, directed=directed); net <- linkEvents(x=net, changeEvent=ev, nodes=nodes)
  dep <- defineDependentEvents(events=ev, nodes=nodes, defaultNetwork=net)
  assign("net", net, globalenv()); assign("dep", dep, globalenv())
  m <- suppressWarnings(estimate(as.formula(paste("dep ~", form)), model="DyNAM", subModel=sub,
         estimationInit=c(list(maxIterations=80, initialDamping=10), init)))
  print(c(coef(m), ll=m$logLikelihood), digits=10) }
t1 <- ev0[,3]; tq <- ceiling(t1/4)
fit(t1, "inertia + recip + indeg", "choice")          # 1.0497675 1.0747602 -0.2707421, ll -38.1300248
fit(t1, "indeg + outdeg", "rate")                     # 0.02063502 0.09236512, ll -55.41977877
fit(t1, "1 + indeg(net) + outdeg(net)", "rate")       # -1.259479036 -0.067194365 0.009146213, ll -94.37558055
fit(tq, "inertia + recip + indeg", "choice")          # 0.9853781955 1.0003945242 -0.3635616050, ll -38.9031366290
fit(tq, "indeg + outdeg", "rate")                     # -0.2559552691 -0.5005636724, ll -54.7027768468
fit(tq, "1 + indeg(net) + outdeg(net)", "rate")       # 0.64072116773 -0.05028715403 -0.27224831861, ll -34.19667944338
fit(t1, "inertia", "choice_coordination", directed=FALSE)   # 0.7290164305, ll -68.8462194594
fit(tq, "inertia", "choice_coordination", directed=FALSE)   # 0.4570586411, ll -70.5803895458
fit(t1, "inertia", "choice")                          # 1.007368923, ll -41.213570853
fit(t1, "inertia(net, window = 5)", "choice")         # 0.92502735, ll -41.15874024 (5.5: same)
fit(t1, "inertia(net, window = 4.5)", "choice")       # 0.8671928, ll -41.76833886
fit(t1, "outdeg(net, window = 5)", "rate")            # -0.1274876, ll -55.34015398
fit(t1, "outdeg(net, window = 4.5)", "rate")          # 0.03669108, ll -55.44260094
fit(t1, "inertia(net, weighted=TRUE)", "choice")      # 0.2889076707, ll -37.7174439566
fit(t1, "indeg(net, weighted=TRUE) + outdeg(net, weighted=TRUE)", "rate")  # -0.08797114084 -0.11872452538, ll -55.01240505209
w <- (1:40 %% 3) + 1
fit(t1, "inertia(net, weighted=TRUE)", "choice", inc=w)   # 0.1389639466, ll -39.2302174436
fit(t1, "indeg(net, weighted=TRUE) + outdeg(net, weighted=TRUE)", "rate", inc=w)  # -0.009244684541 -0.038039028077, ll -55.218661708709
fit(t1, "inertia + recip + indeg + outdeg(net) + trans(net) + cycle(net) + commonSender(net) + commonReceiver(net)", "choice")
  # 2.3222278923 2.6066294426 -1.8832365523 -1.3145112385 1.2728299833 1.1697045224 0.3808007836 -0.6188382957, ll -33.1426170323
opp <- lapply(1:40, function(i) 1:n); opp[[5]] <- 1
fit(t1, "inertia", "choice", init=list(opportunitiesList=opp))   # 1.007368923, ll -40.114958564 (= -41.213570853 + ln 3)
opp <- lapply(1:40, function(i) 1:n); opp[[1]] <- c(1,2,4)
try(fit(t1, "inertia", "choice", init=list(opportunitiesList=opp)))  # ERROR "Active node 3 not available in event 1"
```

The inclusive window rule (an event exactly W ago counts) was confirmed with an own conditional logit in R:
inclusive 5 = goldfish 5, exclusive 5 = 4.5. Hand results and binomial closed forms (see xt_dyn_closedform.do
comments: exp(theta) = a(M-k)/(bk), ll = a ln(a/(Nk)) + b ln(b/(N(M-k))), se = sqrt(N/(ab))) were checked by brute
force in R. Null lls: choice -40 ln 3 = -43.9444915467, rate -40 ln 4 = -55.4517744448.

TODO-R: none open. Every asserted number comes from the R code above or from a hand derivation. The 17-effect and
5-effect rate lists in xt_dyn_options.do check structure only.

## Suspected bugs and code/help or goldfish mismatches

Checked in xt_dyn_suspected.do unless noted otherwise.

1. **S1/S2 (core): large numeric actor ids** get actor index 0. nwset.ado:285-286 (two-mode 132-133) use %9.0g
   labels, while nwfromedge.ado:218, 246 use %18.0g on a float dictionary. This shows for ids > 16,777,216 or
   >= 1e9. It affects nwrem as well.
2. **S3: missing sender or receiver** gives index 0 and a raw Mata error in DynamState. **S4: missing time** is
   accepted, sorted last and used as the last event (with intercept, the waiting time is missing).
3. **S5: opportunities() without the realised receiver is accepted.** goldfish refuses ("Active node 3 not available
   in event 1"). nwdynam.ado:692-701 checks only the ranges.
4. **S6: missing values in the opportunities() variables** pass the range checks (`if !missing()`,
   nwdynam.ado:694, 699). The Mata loop at 704 then subscripts with a missing value.
5. **S7: scalar `__nwdynam_sequential` leaks** after any error raised after nwdynam.ado:657: the 2001 at 659-662, the
   opportunities()/tie()/two-mode checks at 686-772, and the Mata 3498. Only the success path drops it (816).
6. **S8: zero-variance effects are not flagged.** nwdynam.ado:855-860 tests only `missing(se) | se > 100`, and nwrem
   also tests v <= 0. For example, inertiawindow(0.5) on integer times gives a flat likelihood with
   e(unidentified) empty.
7. **S9/S10: sentinel collisions.** inertiawindow(-1) (and the other window options, nwdynam.ado:16, 362-422)
   silently means "no window", while -2 is refused. seed(-1) means "no seed" (642-644).
8. **S11: covariate options are not checked for missing values** (603-640).
9. **S12 (DESIGN, goldfish deviation): weighted effects ignore event values.** goldfish uses them as increments.
   This is reproducible with the R code above (inc = w).
10. **S13 (DESIGN): covariates are read by row position** even when `_nwnode` is present and out of order.
11. **Not tested (needs fault injection):** the `__nwdynam_native` frame is not restored or dropped if optimize()
    errors inside DynamFitUnit1/DynamFitRateUnit1 (unw_dynam.do:457-460 create it, 519-520/683-684 restore it only
    on the normal path). A failing native fit leaves the user in that frame.
12. **Help/code:** nwdynam.sthlp:57 says "most effect sets start the optimizer from random values", but the default
    sets (unit1 engines) use fixed starts and ignore seed() (xt_dyn_estimation.do Block A asserts this).
    e(unidentified) is documented as "missing or above 100", which matches the code but misses item 6.
13. **Effectively unreachable check:** the two-mode sender/receiver direction check (nwdynam.ado:765-772) cannot fire
    through nwset, because nwset twomode always assigns the first variable to mode 1.
14. **Core `_nwsyntax`:** two network names give r(482) "not found" rather than a "one network only" message.
