* dev/dynam_unit21_rateintercept_extras_crosscheck.do -- direct Mata-call
* crosscheck of the WITH-INTERCEPT rate sub-model combined with two-mode
* and weighted=TRUE (batch 18), against real goldfish reference numbers
* captured in dev/dynam_unit21_rateintercept_twomode_crosscheck.R and
* dynam_unit21c_..._weighted_...R.
*
* window() is deliberately NOT tested here (and stays rejected in
* nwdynam.ado): evaluating the fitted log-likelihood directly at real
* goldfish's own reported MLE (not re-optimizing) found a real,
* substantial mismatch (~0.3 in logLik) - traced to goldfish's own
* window() mechanism (createWindowedEvents() in its R source) inserting
* synthetic "dissolve" events into the timeline at event_time + window,
* forcing its piecewise-constant hazard machinery to recompute right at
* each expiry. This engine only recomputes statistics at real dependent
* events, so a windowed contact silently over-counts for the rest of the
* current inter-event interval once intercept's own hazard-INTEGRAL
* aggregation is in play (the no-intercept engine's ordinal partial
* likelihood never integrates over real time, which is why window
* already works correctly there). See dev/dynam_unit21b_
* rateintercept_window_crosscheck.R for the reference numbers that
* exposed this, and nwdynam.ado's own intercept validation block for the
* full account - a real, disclosed architectural gap, not attempted this
* batch.

* (run from the root of this repository)
* The core package nwcommands is found through the environment variable NWCORE (path of the core checkout).
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
* (not needed here: unw_ergm.do belongs to another package)
do unw_dynam.do

di as text "{hline 70}"
di as text "Sub-check 1: intercept + outdeg, TWO-MODE (n1=5, n2=4)"
di as text "{hline 70}"
import delimited "dev/dynam_crosscheck_twomode_intercept_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .)) + 5
mata: events1 = st_data(., ("senderidx","receiveridx","time"))
mata: DynamRateInterceptFitMulti(events1, 9, (0,1,0,0,0), J(1,9,0), J(1,9,0), 1e300, 1e300, 0, 0, 5, "b1", "V1", "ll1")
matrix list b1
scalar d1a = abs(el(b1,1,1) - (-3.52957420))
scalar d1b = abs(el(b1,1,2) - (0.00931169))
scalar d1c = abs(ll1 - (-54.25246354))
di as text "abs diff: intercept=" %9.6f d1a "  outdeg=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001
di as result "  SUB-CHECK 1 PASSED"

di as text "{hline 70}"
di as text "Sub-check 2: intercept + outdeg WEIGHTED=TRUE"
di as text "{hline 70}"
import delimited "dev/dynam_crosscheck_rateintweighted_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata: events3 = st_data(., ("senderidx","receiveridx","time"))
mata: DynamRateInterceptFitMulti(events3, 6, (0,1,0,0,0), J(1,6,0), J(1,6,0), 1e300, 1e300, 0, 1, 0, "b3", "V3", "ll3")
matrix list b3
scalar d3a = abs(el(b3,1,1) - (-2.90509662))
scalar d3b = abs(el(b3,1,2) - (-0.63982127))
scalar d3c = abs(ll3 - (-53.25960504))
di as text "abs diff: intercept=" %9.6f d3a "  outdegweighted=" %9.6f d3b "  logLik=" %9.6f d3c
assert d3a < 0.001
assert d3b < 0.001
assert d3c < 0.001
di as result "  SUB-CHECK 2 PASSED"

di as result "MATCH CONFIRMED (WITH-INTERCEPT rate sub-model reproduces real goldfish exactly combined with two-mode and weighted=TRUE)"
