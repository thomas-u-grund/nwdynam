* dev/dynam_unit3_effectselect_crosscheck.do -- fits unw_dynam.do's own new
* effect-SELECTION engine (the development notes' "effect selection"
* scope item, resolved 2026-09-02 - DynamChoiceFitMulti()/
* DynamRateFitMulti()) against real goldfish subset fits on the SAME
* event data (dev/dynam_crosscheck_events.csv) - matching
* dev/dynam_unit1_choice_crosscheck.do's/dev/dynam_unit2_rate_crosscheck.do's
* own real-data crosscheck discipline (direct Mata calls with n_actors=84
* HARDCODED, matching goldfish's own full 84-actor roster - NOT via
* nwset's own node-count auto-detection, which infers only 54 distinct
* actors from this particular event subset alone and would silently
* fit a materially different, smaller-risk-set model; confirmed directly
* while first drafting this test via nwdynam.ado, which reported "Actors:
* 54" and produced numbers that did not match goldfish at all - not a bug
* in the new effect-selection code, a real actor-count mismatch between
* what nwset can infer from an edgelist alone and goldfish's own
* explicitly-declared 84-actor `actors` data frame), extended here to
* PROPER SUBSETS of each sub-model's own effects, not just the full set.
* nwdynam.ado's own dispatch logic (choosing DynamChoiceFitMulti() for a
* subset vs. the original DynamFitUnit1() for the full set) is certified
* separately, with internally-consistent synthetic data (not needing this
* same 84-vs-54-actor concern), in cscripts/test_nwdynam_ado.do.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439 events -
* see console output captured 2026-09-02):
*   choice ~ inertia + indeg (no recip): inertia=5.960149 (SE 0.173596),
*     indeg=-0.066251 (SE 0.084034), logLik=-733.26
*   choice ~ recip only:                 recip=4.97991 (SE 0.12377), logLik=-1182.8
*   rate ~ outdeg only:                  outdeg=0.435220 (SE 0.020668), logLik=-1783.9

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

di as text "{hline 60}"
di as text "choice ~ inertia + indeg (no recip): active = (1,0,1)"
di as text "{hline 60}"
mata: zero84 = J(1, n_actors, 0)
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero84, zero84, zero84, zero84, zero84, zero84, zero84, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "inertia=" %9.6f el(b1,1,1) "  indeg=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 5.960149)
scalar d2 = abs(el(b1,1,2) - (-0.066251))
scalar d3 = abs(ll1 - (-733.26))
di as text "abs diff: inertia=" %9.6f d1 "  indeg=" %9.6f d2 "  logLik=" %9.4f d3
assert d1 < 0.01
assert d2 < 0.01
assert d3 < 0.05

di as text "{hline 60}"
di as text "choice ~ recip only: active = (0,1,0)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero84, zero84, zero84, zero84, zero84, zero84, zero84, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
assert colsof(b2) == 1
di as text "recip=" %9.6f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d4 = abs(el(b2,1,1) - 4.97991)
scalar d5 = abs(ll2 - (-1182.8))
di as text "abs diff: recip=" %9.6f d4 "  logLik=" %9.4f d5
assert d4 < 0.01
assert d5 < 0.05

di as text "{hline 60}"
di as text "rate ~ outdeg only: active = (0,1)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,1,0,0,0), zero84, zero84, 1e300, 1e300, 0, 0, 0, "b3", "V3", "ll3")
assert colsof(b3) == 1
di as text "outdeg=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d6 = abs(el(b3,1,1) - 0.435220)
scalar d7 = abs(ll3 - (-1783.9))
di as text "abs diff: outdeg=" %9.6f d6 "  logLik=" %9.4f d7
assert d6 < 0.01
assert d7 < 0.05

di as text "{hline 60}"
di as text "Sanity: DynamChoiceFitMulti()/DynamRateFitMulti() with ALL effects active"
di as text "reproduce DynamFitUnit1()/DynamFitRateUnit1()'s own already-verified numbers"
di as text "(both engines climbing the identical likelihood surface when given the same"
di as text "active set, not just individually correct against goldfish)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero84, zero84, zero84, zero84, zero84, zero84, zero84, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
mata: DynamFitUnit1(eventmat, n_actors, "b5", "V5", "ll5")
mata: assert(max(abs(st_matrix("b4") - st_matrix("b5"))) < 0.01)
di as text "  full-set choice: multi=" b4[1,1] " " b4[1,2] " " b4[1,3] "   unit1=" b5[1,1] " " b5[1,2] " " b5[1,3]
di as text "  logLik: multi=" %9.4f ll4 "   unit1=" %9.4f ll5
assert abs(ll4 - ll5) < 0.1

mata: DynamRateFitMulti(eventmat, n_actors, (1,1,0,0,0), zero84, zero84, 1e300, 1e300, 0, 0, 0, "b6", "V6", "ll6")
mata: DynamFitRateUnit1(eventmat, n_actors, "b7", "V7", "ll7")
mata: assert(max(abs(st_matrix("b6") - st_matrix("b7"))) < 0.01)
di as text "  full-set rate:   multi=" b6[1,1] " " b6[1,2] "   unit1=" b7[1,1] " " b7[1,2]
di as text "  logLik: multi=" %9.4f ll6 "   unit1=" %9.4f ll7
assert abs(ll6 - ll7) < 0.1

di as result "MATCH CONFIRMED (subset-effect selection agrees with real goldfish, better than 1e-2 on coefficients, 5e-2 on logLik; full-effect Multi/Unit1 agreement confirmed too)"
