* dev/dynam_unit13_opportunities_crosscheck.do -- opportunitiesList
* support (the development notes' "effect expansion" scope, batch 10,
* resolved 2026-09-02) against real goldfish subset fits on the SAME
* event data - matching dev/dynam_unit12_twomode_crosscheck.do's own
* direct-Mata-call discipline.
*
* MECHANICS (derived from a real toy example first, not guessed - see
* dev/dynam_unit13_opportunities_crosscheck.R's own header comment):
* goldfish's own `opportunitiesList` (choice sub-model only) restricts
* the candidate risk set PER EVENT to a caller-supplied list of
* available node indices. Implemented here as `oppmat', an nevents x n
* 0/1 matrix (row i = which actors are available for event i) added
* on top of the existing per-event `lp' masking (self-exclusion,
* two-mode) - unavailable actors get `lp[j] = lp[j] - 1e300' for that
* event only. `oppmat' with 0 rows (the default, `J(0,0,0)') means "no
* restriction," matching goldfish's own default and reproducing the
* exact pre-existing numbers (verified as a genuine regression check
* below, not just argued equivalent).
*
* Real goldfish reference (8 actors, 30 events - captured 2026-09-02,
* goldfish 1.6.12, see dev/dynam_unit13_opportunities_crosscheck.R for
* the R side):
*   choice ~ inertia + indeg, no opportunitiesList:
*       inertia=0.45766, indeg=-0.54225, logLik=-56.75
*   choice ~ inertia + indeg, WITH per-event opportunitiesList (1
*   random exclusion each event):
*       inertia=0.52986, indeg=-0.52236, logLik=-52.104

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

import delimited "dev/dynam_crosscheck_oppl_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))

mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat = (sender, receiver, eventtime)
n_actors = 8
nevents = rows(eventmat)
zero8 = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "Regression: no opportunitiesList (oppmat = empty) reproduces the exact"
di as text "pre-existing numbers"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
assert colsof(b0) == 2
di as text "inertia=" %9.5f el(b0,1,1) "  indeg=" %9.5f el(b0,1,2) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.45766)
scalar d0b = abs(el(b0,1,2) - (-0.54225))
scalar d0c = abs(ll0 - (-56.75))
di as text "abs diff: inertia=" %9.6f d0a "  indeg=" %9.6f d0b "  logLik=" %9.4f d0c
assert d0a < 0.001
assert d0b < 0.001
assert d0c < 0.01

di as text "{hline 60}"
di as text "choice ~ inertia + indeg, WITH per-event opportunitiesList"
di as text "{hline 60}"

* Build oppmat (nevents x n_actors, 0/1) from the long-format
* availability CSV (one row per (event, available-actor) pair).
preserve
import delimited "dev/dynam_crosscheck_oppl_availability.csv", clear varnames(1)
mata:
oppevent = st_data(., "event")
oppactor = st_data(., "actor")
oppmat = J(nevents, n_actors, 0)
for (oi=1; oi<=rows(oppevent); oi++) {
	oppmat[oppevent[oi], oppactor[oi]] = 1
}
end
restore

mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, oppmat, J(0,0,0), "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "inertia=" %9.5f el(b1,1,1) "  indeg=" %9.5f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - 0.52986)
scalar d1b = abs(el(b1,1,2) - (-0.52236))
scalar d1c = abs(ll1 - (-52.104))
di as text "abs diff: inertia=" %9.6f d1a "  indeg=" %9.6f d1b "  logLik=" %9.4f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.01

di as result "MATCH CONFIRMED (opportunitiesList reproduces real goldfish exactly, both the per-event-restricted fit and the no-restriction regression check)"
