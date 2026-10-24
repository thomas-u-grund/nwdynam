* dev/dynam_unit8_choiceoutdeg_crosscheck.do -- adds the choice
* sub-model's own `outdeg` effect (the development notes' "effect
* expansion" scope, batch 5, resolved 2026-09-02) against real goldfish
* subset fits on the SAME event data - matching
* dev/dynam_unit5_egoalter_crosscheck.do's own direct-Mata-call,
* n_actors=84-hardcoded discipline.
*
* Correcting a real, user-caught documentation error: nwdynam.sthlp's
* own "outdeg: rate sub-model only" label was checked directly against
* goldfish's own doc/goldfishEffects.Rmd and found WRONG - the
* "Structural Effects" table there lists outdeg as valid for BOTH rate
* and choice sub-models. This was a genuine missing effect in the
* choice engine, not a correct scope restriction.
*
* FORMULA - goldfish's own "alter type" outdegree (doc/goldfishEffects.Rmd,
* Outdegree section): s(i,j,t,x) = sum_l I(x(t)_jl > 0), i.e. the
* candidate receiver j's OWN out-degree (how many distinct actors j has
* themselves contacted so far), broadcast with no comparison to the
* sender - mirrors `indeg`'s own "alter type" shape exactly but
* transposed (indeg sums lastcontact's own COLUMN j; outdeg sums its
* own ROW j), reusing the SAME lastcontact matrix already maintained
* for inertia/recip/indeg - no new tracking structure needed.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439
* events - captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit8_choiceoutdeg_crosscheck.R for the R side):
*   choice ~ outdeg:                         outdeg=0.277282, logLik=-1892.9
*   choice ~ inertia + outdeg:               inertia=5.855662, outdeg=0.063952, logLik=-732.37
*   choice ~ inertia + recip + indeg + outdeg: inertia=5.202152, recip=1.920492,
*       indeg=-0.071988, outdeg=-0.254088, logLik=-686.79

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
di as text "choice ~ outdeg (alter type)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
assert colsof(b0) == 1
di as text "outdeg=" %9.5f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.277282)
scalar d0b = abs(ll0 - (-1892.9))
di as text "abs diff: outdeg=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + outdeg"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "inertia=" %9.5f el(b1,1,1) "  outdeg=" %9.5f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 5.855662)
scalar d2 = abs(el(b1,1,2) - 0.063952)
scalar d3 = abs(ll1 - (-732.37))
di as text "abs diff: inertia=" %9.6f d1 "  outdeg=" %9.6f d2 "  logLik=" %9.4f d3
assert d1 < 0.01
assert d2 < 0.01
assert d3 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + recip + indeg + outdeg (all four structural effects together)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
assert colsof(b2) == 4
di as text "inertia=" %9.5f el(b2,1,1) "  recip=" %9.5f el(b2,1,2) "  indeg=" %9.5f el(b2,1,3) "  outdeg=" %9.5f el(b2,1,4) "  logLik=" %9.4f ll2
scalar d4 = abs(el(b2,1,1) - 5.202152)
scalar d5 = abs(el(b2,1,2) - 1.920492)
scalar d6 = abs(el(b2,1,3) - (-0.071988))
scalar d7 = abs(el(b2,1,4) - (-0.254088))
scalar d8 = abs(ll2 - (-686.79))
di as text "abs diff: inertia=" %9.6f d4 "  recip=" %9.6f d5 "  indeg=" %9.6f d6 "  outdeg=" %9.6f d7 "  logLik=" %9.4f d8
assert d4 < 0.01
assert d5 < 0.01
assert d6 < 0.01
assert d7 < 0.01
assert d8 < 0.1

di as result "MATCH CONFIRMED (choice sub-model's own outdeg reproduces real goldfish exactly alone, alongside inertia, and alongside all three other structural effects at once)"
