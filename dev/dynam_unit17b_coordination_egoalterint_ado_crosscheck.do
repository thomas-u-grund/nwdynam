* dev/dynam_unit17b_coordination_egoalterint_ado_crosscheck.do --
* confirms nwdynam.ado's submodel(choice_coordination) egoalterint()
* option (batch 14, 2026-09-02) reproduces dev/dynam_unit17_
* coordination_egoalterint_crosscheck.do's own direct-Mata-call
* numbers through the REAL .ado command, on the SAME data.

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
import delimited "dev/dynam_crosscheck_coord3_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(coordnet3) undirected

nwload coordnet3, xvars
gen floor = 1 in 1
replace floor = 1 in 2
replace floor = 2 in 3
replace floor = 2 in 4
replace floor = 3 in 5
replace floor = 3 in 6
gen floor2 = 5 in 1
replace floor2 = 3 in 2
replace floor2 = 8 in 3
replace floor2 = 1 in 4
replace floor2 = 9 in 5
replace floor2 = 2 in 6

di as text "{hline 60}"
di as text "nwdynam coordnet3, submodel(choice_coordination) inertia egoalterint(floor floor2)"
di as text "{hline 60}"
nwdynam coordnet3, submodel(choice_coordination) inertia egoalterint(floor floor2)
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - (-0.39844340))
scalar d1b = abs(el(b1,1,2) - (-0.03433753))
scalar d1c = abs(ll1 - (-31.60121))
di as text "abs diff: inertia=" %9.6f d1a "  egoAlterInt=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as result "MATCH CONFIRMED (nwdynam.ado's choice_coordination egoalterint() reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
