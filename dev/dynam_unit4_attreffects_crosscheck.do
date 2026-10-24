* dev/dynam_unit4_attreffects_crosscheck.do -- fits unw_dynam.do's own new
* ATTRIBUTE effects (same/diff/sim - the development notes' "effect
* expansion" scope item, first batch, resolved 2026-09-02) against real
* goldfish subset fits on the SAME event data + actor covariate (the
* Social_Evolution dataset's own `floor` variable - which dorm floor each
* actor lives on) - matching dev/dynam_unit3_effectselect_crosscheck.do's
* own direct-Mata-call, n_actors=84-hardcoded discipline (see that file's
* own comment for why nwset's own node-count auto-detection cannot be
* used here either).
*
* Real goldfish reference (Social_Evolution$calls + actors$floor, 84
* actors, 439 events - see console output captured 2026-09-02):
*   choice ~ same(floor) alone: same=0.67677 (SE 0.11748), logLik=-1925.4
*   choice ~ diff(floor) alone: diff=-0.529022 (SE 0.035172), logLik=-1775.6
*   choice ~ sim(floor) alone:  sim=0.529022 (SE 0.035172), logLik=-1775.6
*     (sim = -diff exactly, by construction - goldfish's own documented
*     relationship, not a coincidence - see unw_dynam.do's own comment)
*   choice ~ inertia + same(floor): inertia=5.95911 (SE 0.15583),
*     same=-0.37779 (SE 0.15022), logLik=-730.3

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

mata: zero6 = J(1,6,0)

di as text "{hline 60}"
di as text "choice ~ same(floor) alone: active = (0,0,0,1,0,0)"
di as text "{hline 60}"
mata: zero84b = J(1, n_actors, 0)
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), floorvec, floorvec, floorvec, zero84b, zero84b, zero84b, zero84b, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "same=" %9.6f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 0.67677)
scalar d2 = abs(ll1 - (-1925.4))
di as text "abs diff: same=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.01
assert d2 < 0.1

di as text "{hline 60}"
di as text "choice ~ diff(floor) alone: active = (0,0,0,0,1,0)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0), floorvec, floorvec, floorvec, zero84b, zero84b, zero84b, zero84b, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
di as text "diff=" %9.6f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - (-0.529022))
scalar d4 = abs(ll2 - (-1775.6))
di as text "abs diff: diff=" %9.6f d3 "  logLik=" %9.4f d4
assert d3 < 0.01
assert d4 < 0.1

di as text "{hline 60}"
di as text "choice ~ sim(floor) alone: active = (0,0,0,0,0,1)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0), floorvec, floorvec, floorvec, zero84b, zero84b, zero84b, zero84b, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "sim=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d5 = abs(el(b3,1,1) - 0.529022)
scalar d6 = abs(ll3 - (-1775.6))
di as text "abs diff: sim=" %9.6f d5 "  logLik=" %9.4f d6
assert d5 < 0.01
assert d6 < 0.1

* sim = -diff exactly, on THIS engine too, not just goldfish's own
* documented claim about itself - a free internal identity check. Not
* bit-exact: sim and diff are each fit by an INDEPENDENT BFGS run (from
* random starting points, on a genuinely reparameterized - theta ->
* -theta - likelihood surface), so a tiny optimizer-path difference
* between the two runs is expected, not a formula bug; same tolerance
* class as the goldfish-agreement checks above, not machine epsilon.
scalar d7 = abs(el(b3,1,1) - (-el(b2,1,1)))
di as text "internal check: sim vs -diff, abs diff=" %9.8f d7
assert d7 < 0.001
assert abs(ll3 - ll2) < 0.001

di as text "{hline 60}"
di as text "choice ~ inertia + same(floor): active = (1,0,0,1,0,0)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), floorvec, floorvec, floorvec, zero84b, zero84b, zero84b, zero84b, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
assert colsof(b4) == 2
di as text "inertia=" %9.5f el(b4,1,1) "  same=" %9.5f el(b4,1,2) "  logLik=" %9.4f ll4
scalar d8 = abs(el(b4,1,1) - 5.95911)
scalar d9 = abs(el(b4,1,2) - (-0.37779))
scalar d10 = abs(ll4 - (-730.3))
di as text "abs diff: inertia=" %9.6f d8 "  same=" %9.6f d9 "  logLik=" %9.4f d10
assert d8 < 0.01
assert d9 < 0.01
assert d10 < 0.1

di as result "MATCH CONFIRMED (same/diff/sim attribute effects agree with real goldfish, better than 1e-2 on coefficients, 0.1 on logLik; sim=-diff identity confirmed internally)"
