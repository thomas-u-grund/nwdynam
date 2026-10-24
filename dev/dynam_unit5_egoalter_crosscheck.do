* dev/dynam_unit5_egoalter_crosscheck.do -- fits unw_dynam.do's own new
* ego()/alter() attribute effects (the development notes' "effect
* expansion" scope item, batch 2, resolved 2026-09-02) against real
* goldfish subset fits on the SAME event data + actor covariate (the
* Social_Evolution dataset's own `floor` variable) - matching
* dev/dynam_unit4_attreffects_crosscheck.do's own direct-Mata-call,
* n_actors=84-hardcoded discipline.
*
* ego(z) = z_i (rate sub-model only - the candidate SENDER's own static
* covariate, no comparison to anyone); alter(z) = z_j (choice sub-model
* only - the candidate RECEIVER's own static covariate) - the simplest
* possible attribute effects in goldfish's own catalog, formulas read
* directly from doc/goldfishEffects.Rmd.
*
* Real goldfish reference (Social_Evolution$calls + actors$floor, 84
* actors, 439 events - see console output captured 2026-09-02):
*   rate ~ ego(floor) alone:            ego=-0.351453 (SE 0.024845), logLik=-1817
*   rate ~ ego(floor) + outdeg:         ego=-0.399894, outdeg=0.420608, logLik=-1669.9
*   choice ~ alter(floor) alone:        alter=-0.304892 (SE 0.023538), logLik=-1837.4
*   choice ~ inertia + alter(floor):    inertia=5.686052, alter=-0.166814, logLik=-721.17

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
end

preserve
import delimited "dev/dynam_crosscheck_actorcov.csv", clear varnames(1)
sort actor
assert _N == 84
mata: floorvec = st_data(1::84, "floor")'
restore

mata: zero7 = J(1, n_actors, 0)

di as text "{hline 60}"
di as text "rate ~ ego(floor) alone: active = (0,0,1)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,0,1,0,0), floorvec, floorvec, 1e300, 1e300, 0, 0, 0, "b1", "V1", "ll1")
di as text "ego=" %9.6f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - (-0.351453))
scalar d2 = abs(ll1 - (-1817))
di as text "abs diff: ego=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.01
assert d2 < 0.1

di as text "{hline 60}"
di as text "rate ~ ego(floor) + outdeg: active = (0,1,1)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,1,1,0,0), floorvec, floorvec, 1e300, 1e300, 0, 0, 0, "b2", "V2", "ll2")
assert colsof(b2) == 2
di as text "outdeg=" %9.6f el(b2,1,1) "  ego=" %9.6f el(b2,1,2) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - 0.420608)
scalar d4 = abs(el(b2,1,2) - (-0.399894))
scalar d5 = abs(ll2 - (-1669.9))
di as text "abs diff: outdeg=" %9.6f d3 "  ego=" %9.6f d4 "  logLik=" %9.4f d5
assert d3 < 0.01
assert d4 < 0.01
assert d5 < 0.1

di as text "{hline 60}"
di as text "choice ~ alter(floor) alone: active = (0,0,0,0,0,0,1)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0), zero7, zero7, zero7, floorvec, zero7, zero7, zero7, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "alter=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d6 = abs(el(b3,1,1) - (-0.304892))
scalar d7 = abs(ll3 - (-1837.4))
di as text "abs diff: alter=" %9.6f d6 "  logLik=" %9.4f d7
assert d6 < 0.01
assert d7 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + alter(floor): active = (1,0,0,0,0,0,1)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0), zero7, zero7, zero7, floorvec, zero7, zero7, zero7, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
assert colsof(b4) == 2
di as text "inertia=" %9.5f el(b4,1,1) "  alter=" %9.5f el(b4,1,2) "  logLik=" %9.4f ll4
scalar d8 = abs(el(b4,1,1) - 5.686052)
scalar d9 = abs(el(b4,1,2) - (-0.166814))
scalar d10 = abs(ll4 - (-721.17))
di as text "abs diff: inertia=" %9.6f d8 "  alter=" %9.6f d9 "  logLik=" %9.4f d10
assert d8 < 0.01
assert d9 < 0.01
assert d10 < 0.1

di as result "MATCH CONFIRMED (ego/alter attribute effects agree with real goldfish, better than 1e-2 on coefficients, 0.1 on logLik)"
