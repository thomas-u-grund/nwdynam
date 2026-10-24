* dev/dynam_unit19_coordination_tie_crosscheck.do --
* choice_coordination effect expansion, batch 16, resolved 2026-09-02 -
* extends tie() (cross-network effects, v1 scope) to
* choice_coordination, verified against real goldfish subset fits on
* an undirected 6-actor/12-event toy network with a separate static
* exogenous network.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit19_coordination_tie_crosscheck.R for the R side):
*   tie(staticnet) alone:            tie=0.7081185, logLik=-31.79857
*   inertia + tie(staticnet):        inertia=-0.3867456, tie=0.6756129,
*                                     logLik=-31.43784

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

import delimited "dev/dynam_crosscheck_coordtie_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat = (sender, receiver, eventtime)
n_actors = 6
zero6 = J(1, n_actors, 0)
end

preserve
import delimited "dev/dynam_crosscheck_coordtie_staticmat.csv", clear varnames(1)
mata: tiemat = st_data(., .)
restore

di as text "{hline 60}"
di as text "choice_coordination ~ tie(staticnet)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,1), zero6, zero6, zero6, zero6, zero6, zero6, zero6, tiemat, "b0", "V0", "ll0")
di as text "tie=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.7081185)
scalar d0b = abs(ll0 - (-31.79857))
di as text "abs diff: tie=" %9.6f d0a "  logLik=" %9.6f d0b
assert d0a < 0.001
assert d0b < 0.001

di as text "{hline 60}"
di as text "choice_coordination ~ inertia + tie(staticnet)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,1), zero6, zero6, zero6, zero6, zero6, zero6, zero6, tiemat, "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "inertia=" %9.6f el(b1,1,1) "  tie=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - (-0.3867456))
scalar d1b = abs(el(b1,1,2) - 0.6756129)
scalar d1c = abs(ll1 - (-31.43784))
di as text "abs diff: inertia=" %9.6f d1a "  tie=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as result "MATCH CONFIRMED (choice_coordination's tie() reproduces real goldfish exactly, alone and combined with inertia)"
