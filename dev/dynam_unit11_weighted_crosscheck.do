* dev/dynam_unit11_weighted_crosscheck.do -- adds weighted=TRUE support
* for inertia/recip/indeg/outdeg (the development notes' "effect
* expansion" scope, batch 8, resolved 2026-09-02, continuing "work
* your way through the list") against real goldfish subset fits on the
* SAME event data - matching dev/dynam_unit10_tertius_egoalterint_crosscheck.do's
* own direct-Mata-call, n_actors=84-hardcoded discipline.
*
* SEMANTICS - goldfish's own `weighted=TRUE` argument: the effect
* reads the cumulative COUNT of prior events for that dyad/actor
* (repeated calls between the same sender/receiver all count) instead
* of binary presence (ever contacted, yes/no). Verified only for
* inertia/recip/indeg/outdeg in this batch - goldfish's own effect
* table also lists nodeTrans/tertius/tie as weighted-capable, but
* those remain a disclosed follow-on (tie is not implemented at all
* yet; nodeTrans/tertius would each need their own weighted matrix
* products, not yet derived/verified). Mutually exclusive with that
* same effect's own window in this implementation - not verified
* together, and nwdynam.ado rejects the combination outright.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439
* events - captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit11_weighted_crosscheck.R for the R side):
*   choice ~ inertia(weighted=TRUE):              inertia=0.451967, logLik=-1130.5
*   choice ~ recip(weighted=TRUE):                recip=0.294334,   logLik=-1363.1
*   choice ~ indeg(weighted=TRUE):                indeg=0.0473282,  logLik=-1649.0
*   choice ~ outdeg(weighted=TRUE):                outdeg=0.0438901, logLik=-1704.4
*   rate ~ indeg(weighted=TRUE):                  indeg=0.0410526,  logLik=-1724.4
*   rate ~ outdeg(weighted=TRUE):                 outdeg=0.044179,  logLik=-1651.2
*   choice ~ inertia(weighted=TRUE)+recip(weighted=TRUE):
*       inertia=0.396725, recip=0.062008, logLik=-1120.0

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
zero8 = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "choice ~ inertia(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 1, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
di as text "inertia=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.451967)
scalar d0b = abs(ll0 - (-1130.5))
di as text "abs diff: inertia=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "choice ~ recip(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 1, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "recip=" %9.6f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 0.294334)
scalar d2 = abs(ll1 - (-1363.1))
di as text "abs diff: recip=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.01
assert d2 < 0.1

di as text "{hline 60}"
di as text "choice ~ indeg(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 1, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
di as text "indeg=" %9.6f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - 0.0473282)
scalar d4 = abs(ll2 - (-1649.0))
di as text "abs diff: indeg=" %9.6f d3 "  logLik=" %9.4f d4
assert d3 < 0.01
assert d4 < 0.1

di as text "{hline 60}"
di as text "choice ~ outdeg(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 1, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "outdeg=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d5 = abs(el(b3,1,1) - 0.0438901)
scalar d6 = abs(ll3 - (-1704.4))
di as text "abs diff: outdeg=" %9.6f d5 "  logLik=" %9.4f d6
assert d5 < 0.01
assert d6 < 0.1

di as text "{hline 60}"
di as text "rate ~ indeg(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (1,0,0,0,0), zero8, zero8, 1e300, 1e300, 1, 0, 0, "b4", "V4", "ll4")
di as text "indeg=" %9.6f el(b4,1,1) "  logLik=" %9.4f ll4
scalar d7 = abs(el(b4,1,1) - 0.0410526)
scalar d8 = abs(ll4 - (-1724.4))
di as text "abs diff: indeg=" %9.6f d7 "  logLik=" %9.4f d8
assert d7 < 0.01
assert d8 < 0.1

di as text "{hline 60}"
di as text "rate ~ outdeg(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,1,0,0,0), zero8, zero8, 1e300, 1e300, 0, 1, 0, "b5", "V5", "ll5")
di as text "outdeg=" %9.6f el(b5,1,1) "  logLik=" %9.4f ll5
scalar d9 = abs(el(b5,1,1) - 0.044179)
scalar d10 = abs(ll5 - (-1651.2))
di as text "abs diff: outdeg=" %9.6f d9 "  logLik=" %9.4f d10
assert d9 < 0.01
assert d10 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia(weighted=TRUE) + recip(weighted=TRUE)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 1, 1, 0, 0, 0, J(0,0,0), J(0,0,0), "b6", "V6", "ll6")
assert colsof(b6) == 2
di as text "inertia=" %9.6f el(b6,1,1) "  recip=" %9.6f el(b6,1,2) "  logLik=" %9.4f ll6
scalar d11 = abs(el(b6,1,1) - 0.396725)
scalar d12 = abs(el(b6,1,2) - 0.062008)
scalar d13 = abs(ll6 - (-1120.0))
di as text "abs diff: inertia=" %9.6f d11 "  recip=" %9.6f d12 "  logLik=" %9.4f d13
assert d11 < 0.01
assert d12 < 0.01
assert d13 < 0.1

di as text "{hline 60}"
di as text "Regression: unweighted inertia unaffected by weighted-flag plumbing"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b7", "V7", "ll7")
di as text "inertia=" %9.4f el(b7,1,1) "  logLik=" %9.4f ll7
scalar d14 = abs(el(b7,1,1) - 5.8971)
scalar d15 = abs(ll7 - (-733.58))
di as text "abs diff: inertia=" %9.6f d14 "  logLik=" %9.4f d15
assert d14 < 0.01
assert d15 < 0.1

di as result "MATCH CONFIRMED (weighted=TRUE reproduces real goldfish for inertia/recip/indeg/outdeg under both sub-models where applicable, alone and combined; unweighted case unaffected)"
