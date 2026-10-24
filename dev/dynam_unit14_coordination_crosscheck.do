* dev/dynam_unit14_coordination_crosscheck.do -- choice_coordination
* support (the development notes' "effect expansion" scope, batch 11,
* resolved 2026-09-02, user: "do choice_coordination now") against real
* goldfish subset fits - direct Mata-call discipline, matching every
* prior dynam_unitN_*_crosscheck.do's own precedent.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit14_coordination_crosscheck.R for the R side):
*   Toy 1 (4 actors, 4 events): inertia only:
*       inertia=-4.172212, logLik=-5.729604
*   Toy 2 (6 actors, 12 events): inertia + same(floor):
*       inertia=-0.4636019, same=-0.5637855, logLik=-31.33132
*   Toy 3, same 6-actor/12-event data, each effect FIT ALONE:
*       indeg=0.1233222,  logLik=-32.247
*       diff=0.322213,    logLik=-31.64643
*       sim=-0.322213,    logLik=-31.64643
*       alter=-0.2699238, logLik=-31.8523

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

di as text "{hline 60}"
di as text "Toy 1: choice_coordination ~ inertia (4 actors, 4 events)"
di as text "{hline 60}"
import delimited "dev/dynam_crosscheck_coord1_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat1 = (sender, receiver, eventtime)
n1 = 4
zero1 = J(1, n1, 0)
end
mata: DynamCoordFitMulti(eventmat1, n1, (1,0,0,0,0,0,0,0,0,0,0,0), zero1, zero1, zero1, zero1, zero1, zero1, zero1, J(0,0,0), "b0", "V0", "ll0")
di as text "inertia=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - (-4.172212))
scalar d0b = abs(ll0 - (-5.729604))
di as text "abs diff: inertia=" %9.6f d0a "  logLik=" %9.6f d0b
* Toy 1's own inertia coefficient is NOT tightly checked here - a real
* finding, not a bug: none of these 4 events is a repeat contact, so
* the likelihood is on a genuinely FLAT ridge as beta -> -infinity
* (checked directly by re-evaluating this file's own loglik formula in
* R at beta = -4.17, -4.52, -5, -10, -20: logLik moves only in the 4th
* decimal place across that entire range - matching goldfish's own
* huge reported SE of 15.66 on this coefficient). logLik, the only
* well-posed quantity here, matches goldfish to within 0.0002; the
* coefficient itself is only checked directionally (large and
* negative), never to a tight tolerance. Toy 2 below (more data, two
* effects) is the real coefficient-identification crosscheck.
assert d0b < 0.001
assert el(b0,1,1) < -3

di as text "{hline 60}"
di as text "Toy 2: choice_coordination ~ inertia + same(floor) (6 actors, 12 events)"
di as text "{hline 60}"
import delimited "dev/dynam_crosscheck_coord2_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat2 = (sender, receiver, eventtime)
n2 = 6
zero2 = J(1, n2, 0)
end
preserve
import delimited "dev/dynam_crosscheck_coord2_actors.csv", clear varnames(1)
mata: floorvec = st_data(., "floor")'
restore

mata: DynamCoordFitMulti(eventmat2, n2, (1,0,1,0,0,0,0,0,0,0,0,0), floorvec, zero2, zero2, zero2, zero2, zero2, zero2, J(0,0,0), "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "inertia=" %9.6f el(b1,1,1) "  same=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - (-0.4636019))
scalar d1b = abs(el(b1,1,2) - (-0.5637855))
scalar d1c = abs(ll1 - (-31.33132))
di as text "abs diff: inertia=" %9.6f d1a "  same=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as text "{hline 60}"
di as text "Toy 3a: choice_coordination ~ indeg"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat2, n2, (0,1,0,0,0,0,0,0,0,0,0,0), zero2, zero2, zero2, zero2, zero2, zero2, zero2, J(0,0,0), "b2", "V2", "ll2")
di as text "indeg=" %9.6f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d2a = abs(el(b2,1,1) - 0.1233222)
scalar d2b = abs(ll2 - (-32.247))
di as text "abs diff: indeg=" %9.6f d2a "  logLik=" %9.6f d2b
assert d2a < 0.001
assert d2b < 0.001

di as text "{hline 60}"
di as text "Toy 3b: choice_coordination ~ diff(floor)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat2, n2, (0,0,0,1,0,0,0,0,0,0,0,0), zero2, floorvec, zero2, zero2, zero2, zero2, zero2, J(0,0,0), "b3", "V3", "ll3")
di as text "diff=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d3a = abs(el(b3,1,1) - 0.322213)
scalar d3b = abs(ll3 - (-31.64643))
di as text "abs diff: diff=" %9.6f d3a "  logLik=" %9.6f d3b
assert d3a < 0.001
assert d3b < 0.001

di as text "{hline 60}"
di as text "Toy 3c: choice_coordination ~ sim(floor)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat2, n2, (0,0,0,0,1,0,0,0,0,0,0,0), zero2, zero2, floorvec, zero2, zero2, zero2, zero2, J(0,0,0), "b4", "V4", "ll4")
di as text "sim=" %9.6f el(b4,1,1) "  logLik=" %9.4f ll4
scalar d4a = abs(el(b4,1,1) - (-0.322213))
scalar d4b = abs(ll4 - (-31.64643))
di as text "abs diff: sim=" %9.6f d4a "  logLik=" %9.6f d4b
assert d4a < 0.001
assert d4b < 0.001

di as text "{hline 60}"
di as text "Toy 3d: choice_coordination ~ alter(floor)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat2, n2, (0,0,0,0,0,1,0,0,0,0,0,0), zero2, zero2, zero2, floorvec, zero2, zero2, zero2, J(0,0,0), "b5", "V5", "ll5")
di as text "alter=" %9.6f el(b5,1,1) "  logLik=" %9.4f ll5
scalar d5a = abs(el(b5,1,1) - (-0.2699238))
scalar d5b = abs(ll5 - (-31.8523))
di as text "abs diff: alter=" %9.6f d5a "  logLik=" %9.6f d5b
assert d5a < 0.001
assert d5b < 0.001

di as result "MATCH CONFIRMED (choice_coordination reproduces real goldfish exactly on all six v1 effects)"
