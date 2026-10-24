* dev/dynam_unit12_twomode_crosscheck.do -- two-mode (bipartite) DyNAM
* support (the development notes' "effect expansion" scope, batch 9,
* resolved 2026-09-02, user: "continue") against real goldfish subset
* fits on the SAME event data - matching
* dev/dynam_unit11_weighted_crosscheck.do's own direct-Mata-call
* discipline, but with a genuinely different n1/n2/eventmat shape
* (6 people + 4 orgs = 10 actors total, not the usual 84-actor
* Social_Evolution dataset).
*
* See unw_dynam.do's own header comment (above dynam_choice_loglik_
* grad_multi()) for the full two-mode investigation: which effects
* goldfish itself hard-rejects for two-mode, and the deeper finding
* that even simple attribute effects (same/ego) compute something
* goldfish-internal-specific that could not be confidently
* reverse-engineered - so only inertia/indeg(choice)/outdeg(rate) are
* verified and shipped here.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit12_twomode_crosscheck.R for the R side):
*   choice ~ inertia:            inertia=-0.42729, logLik=-82.591
*   choice ~ indeg:               indeg=0.12673,   logLik=-82.692
*   choice ~ inertia + indeg:    inertia=-0.73342, indeg=0.23018, logLik=-81.324
*   rate ~ outdeg:                outdeg=0.23797,  logLik=-106.63

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

import delimited "dev/dynam_crosscheck_twomode_events.csv", clear varnames(1)

* sender/receiver are string labels ("P1".."P6", "O1".."O4") - build
* the combined 1..10 actor index directly (mode 1 = people = 1..6,
* mode 2 = orgs = 7..10), matching nwset's own contiguous convention
* confirmed directly against real nwset output (see
* cscripts/test_nwdynam_twomode_ado.do's own header comment).
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .)) + 6

mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat = (sender, receiver, eventtime)
n_actors = 10
n1 = 6
zero10 = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "choice ~ inertia (two-mode, n1=6)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero10, zero10, zero10, zero10, zero10, zero10, zero10, 1e300, 1e300, 0, 0, 0, 0, n1, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
di as text "inertia=" %9.5f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - (-0.42729))
scalar d0b = abs(ll0 - (-82.591))
di as text "abs diff: inertia=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.001
assert d0b < 0.01

di as text "{hline 60}"
di as text "choice ~ indeg (two-mode, n1=6)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero10, zero10, zero10, zero10, zero10, zero10, zero10, 1e300, 1e300, 0, 0, 0, 0, n1, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "indeg=" %9.5f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 0.12673)
scalar d2 = abs(ll1 - (-82.692))
di as text "abs diff: indeg=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.001
assert d2 < 0.01

di as text "{hline 60}"
di as text "choice ~ inertia + indeg (two-mode, n1=6)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero10, zero10, zero10, zero10, zero10, zero10, zero10, 1e300, 1e300, 0, 0, 0, 0, n1, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
assert colsof(b2) == 2
di as text "inertia=" %9.5f el(b2,1,1) "  indeg=" %9.5f el(b2,1,2) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - (-0.73342))
scalar d4 = abs(el(b2,1,2) - 0.23018)
scalar d5 = abs(ll2 - (-81.324))
di as text "abs diff: inertia=" %9.6f d3 "  indeg=" %9.6f d4 "  logLik=" %9.4f d5
assert d3 < 0.001
assert d4 < 0.001
assert d5 < 0.01

di as text "{hline 60}"
di as text "rate ~ outdeg (two-mode, n1=6)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,1,0,0,0), zero10, zero10, 1e300, 1e300, 0, 0, n1, "b3", "V3", "ll3")
di as text "outdeg=" %9.5f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d6 = abs(el(b3,1,1) - 0.23797)
scalar d7 = abs(ll3 - (-106.63))
di as text "abs diff: outdeg=" %9.6f d6 "  logLik=" %9.4f d7
assert d6 < 0.001
assert d7 < 0.01

di as text "{hline 60}"
di as text "Regression: n1=0 (one-mode) path unaffected by the two-mode masking code"
di as text "{hline 60}"
import delimited "dev/dynam_crosscheck_events.csv", clear varnames(1)
mata:
sender84 = st_data(., "sender")
receiver84 = st_data(., "receiver")
eventtime84 = st_data(., "eventtime")
eventmat84 = (sender84, receiver84, eventtime84)
n84 = 84
zero84c = J(1, n84, 0)
end
mata: DynamChoiceFitMulti(eventmat84, n84, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero84c, zero84c, zero84c, zero84c, zero84c, zero84c, zero84c, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
di as text "inertia=" %9.4f el(b4,1,1) "  logLik=" %9.4f ll4
scalar d8 = abs(el(b4,1,1) - 5.8971)
scalar d9 = abs(ll4 - (-733.58))
di as text "abs diff: inertia=" %9.6f d8 "  logLik=" %9.4f d9
assert d8 < 0.01
assert d9 < 0.1

di as result "MATCH CONFIRMED (two-mode inertia/indeg-choice/outdeg-rate reproduce real goldfish exactly; one-mode path unaffected)"
