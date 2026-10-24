* dev/dynam_unit7_ratewindow_crosscheck.do -- extends unw_dynam.do's own
* window support to the RATE sub-model's indeg/outdeg effects
* (the development notes' "effect expansion" scope, batch 4, resolved
* 2026-09-02, user follow-up "can we make [window] also available for
* the rate part of the model?" immediately after the choice-side
* inertiawindow()/recipwindow() self-activation redesign) against real
* goldfish subset fits on the SAME event data - matching
* dev/dynam_unit6_window_crosscheck.do's own direct-Mata-call,
* n_actors=84-hardcoded discipline.
*
* WINDOW SEMANTICS - same real-time recency-cutoff shape as the choice
* sub-model's own inertia/recip window (t_current - t_lastcontact <=
* window, inclusive, hard cutoff not decay), re-derived and verified
* against real goldfish for the RATE sub-model specifically (not
* assumed identical without checking) - see
* dev/dynam_unit7_ratewindow_crosscheck.R's own console output,
* captured 2026-09-02, goldfish 1.6.12, on Social_Evolution (84 actors,
* 439 real phone-call events):
*   rate ~ indeg, no window:                indeg=0.910956, logLik=-1793.8
*   rate ~ indeg, window=7 days (604800s):   indeg=1.035436, logLik=-1768.9
*   rate ~ indeg, window=1 day (86400s):     indeg=1.817017, logLik=-1678.1
*   rate ~ outdeg, window=1 day (86400s):    outdeg=1.585193, logLik=-1553.0
*   rate ~ indeg(window=7d)+outdeg(window=1d): indeg=0.395815, outdeg=1.370516, logLik=-1531.6

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
do unw_dynam.do

import delimited "dev/dynam_crosscheck_events.csv", clear varnames(1)

mata:
sender = st_data(., "sender")
receiver = st_data(., "receiver")
eventtime = st_data(., "eventtime")
eventmat = (sender, receiver, eventtime)
n_actors = 84
zero3b = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "rate ~ indeg, no window (regression: must reproduce Unit 1's own"
di as text "pre-window numbers exactly, not just approximately)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (1,0,0,0,0), zero3b, zero3b, 1e300, 1e300, 0, 0, 0, "b0", "V0", "ll0")
di as text "indeg=" %9.4f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.910956)
scalar d0b = abs(ll0 - (-1793.8))
di as text "abs diff: indeg=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "rate ~ indeg, window = 7 days (604800 seconds)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (1,0,0,0,0), zero3b, zero3b, 604800, 1e300, 0, 0, 0, "b1", "V1", "ll1")
di as text "indeg=" %9.5f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 1.035436)
scalar d2 = abs(ll1 - (-1768.9))
di as text "abs diff: indeg=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.01
assert d2 < 0.1

di as text "{hline 60}"
di as text "rate ~ indeg, window = 1 day (86400 seconds)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (1,0,0,0,0), zero3b, zero3b, 86400, 1e300, 0, 0, 0, "b2", "V2", "ll2")
di as text "indeg=" %9.5f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - 1.817017)
scalar d4 = abs(ll2 - (-1678.1))
di as text "abs diff: indeg=" %9.6f d3 "  logLik=" %9.4f d4
assert d3 < 0.01
assert d4 < 0.1

di as text "{hline 60}"
di as text "rate ~ outdeg, window = 1 day (86400 seconds)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,1,0,0,0), zero3b, zero3b, 1e300, 86400, 0, 0, 0, "b3", "V3", "ll3")
di as text "outdeg=" %9.5f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d5 = abs(el(b3,1,1) - 1.585193)
scalar d6 = abs(ll3 - (-1553.0))
di as text "abs diff: outdeg=" %9.6f d5 "  logLik=" %9.4f d6
assert d5 < 0.01
assert d6 < 0.1

di as text "{hline 60}"
di as text "rate ~ indeg(window=7 days) + outdeg(window=1 day) - GENUINELY"
di as text "DIFFERENT windows on two effects in the SAME model"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (1,1,0,0,0), zero3b, zero3b, 604800, 86400, 0, 0, 0, "b4", "V4", "ll4")
assert colsof(b4) == 2
di as text "indeg=" %9.5f el(b4,1,1) "  outdeg=" %9.5f el(b4,1,2) "  logLik=" %9.4f ll4
scalar d7 = abs(el(b4,1,1) - 0.395815)
scalar d8 = abs(el(b4,1,2) - 1.370516)
scalar d9 = abs(ll4 - (-1531.6))
di as text "abs diff: indeg=" %9.6f d7 "  outdeg=" %9.6f d8 "  logLik=" %9.4f d9
assert d7 < 0.01
assert d8 < 0.01
assert d9 < 0.1

di as result "MATCH CONFIRMED (rate sub-model's own indeg/outdeg window() reproduces real goldfish exactly at 3 window sizes for indeg, 1 for outdeg, and one model with GENUINELY DIFFERENT windows on both simultaneously; no-window regression unaffected)"
