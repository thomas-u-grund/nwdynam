* dev/dynam_unit14b_coordination_ado_crosscheck.do -- confirms the
* nwdynam.ado submodel(choice_coordination) option (batch 11,
* 2026-09-02) reproduces dev/dynam_unit14_coordination_crosscheck.do's
* own direct-Mata-call numbers when run through the REAL .ado command,
* on the SAME data (dev/dynam_crosscheck_coord2_events.csv, built by
* dev/dynam_unit14_coordination_crosscheck.R against real goldfish
* 1.6.12 - actors 1-6 have floor = 1,1,2,2,3,3, matching that R
* script's own `actors2' data frame).

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
* (not needed here: unw_ergm.do belongs to another package)
do unw_dynam.do
do nwdynam.ado

nwclear
import delimited "dev/dynam_crosscheck_coord2_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(coordnet) undirected

di as text "{hline 60}"
di as text "Test: submodel(choice_coordination) requires undirected network"
di as text "{hline 60}"
import delimited "dev/dynam_crosscheck_coord2_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(directedver)
capture noisily nwdynam directedver, submodel(choice_coordination)
assert _rc == 198
di as text "  correctly rejected directed network (rc=" _rc ")"

di as text "{hline 60}"
di as text "nwdynam coordnet, submodel(choice_coordination) inertia same(floor)"
di as text "{hline 60}"
nwload coordnet, xvars
gen floor = 1 in 1
replace floor = 1 in 2
replace floor = 2 in 3
replace floor = 2 in 4
replace floor = 3 in 5
replace floor = 3 in 6

nwdynam coordnet, submodel(choice_coordination) inertia same(floor)
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - (-0.4636019))
scalar d1b = abs(el(b1,1,2) - (-0.5637855))
scalar d1c = abs(ll1 - (-31.33132))
di as text "abs diff: inertia=" %9.6f d1a "  same=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as text "{hline 60}"
di as text "Test: choice-only/rate-only effects rejected under choice_coordination"
di as text "{hline 60}"
capture noisily nwdynam coordnet, submodel(choice_coordination) recip
assert _rc == 198
di as text "  correctly rejected recip (rc=" _rc ")"
capture noisily nwdynam coordnet, submodel(choice_coordination) cycle
assert _rc == 198
di as text "  correctly rejected cycle (rc=" _rc ")"
capture noisily nwdynam coordnet, submodel(choice_coordination) ego(floor)
assert _rc == 198
di as text "  correctly rejected ego() (rc=" _rc ")"

di as text "{hline 60}"
di as text "Test: opportunities() rejected under choice_coordination"
di as text "{hline 60}"
capture noisily nwdynam coordnet, submodel(choice_coordination) opportunities(floor floor)
assert _rc == 198
di as text "  correctly rejected opportunities() (rc=" _rc ")"

di as result "MATCH CONFIRMED (nwdynam.ado's submodel(choice_coordination) reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
