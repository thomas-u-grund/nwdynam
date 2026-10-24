* dev/dynam_unit16b_coordination_tertius_four_ado_crosscheck.do --
* confirms nwdynam.ado's submodel(choice_coordination) tertius()/four
* options (batch 13, 2026-09-02) reproduce dev/dynam_unit16_
* coordination_tertius_four_crosscheck.do's own direct-Mata-call
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

di as text "{hline 60}"
di as text "nwdynam coordnet3, submodel(choice_coordination) inertia tertius(floor) four"
di as text "{hline 60}"
nwdynam coordnet3, submodel(choice_coordination) inertia tertius(floor) four
matrix b2 = e(b)
scalar ll2 = e(ll)
scalar d2a = abs(el(b2,1,1) - (-0.57080704))
scalar d2b = abs(el(b2,1,2) - 0.32457616)
scalar d2c = abs(el(b2,1,3) - 0.01745949)
scalar d2d = abs(ll2 - (-31.10879))
di as text "abs diff: inertia=" %9.6f d2a "  tertius=" %9.6f d2b "  four=" %9.6f d2c "  logLik=" %9.6f d2d
assert d2a < 0.001
assert d2b < 0.001
assert d2c < 0.001
assert d2d < 0.001

di as result "MATCH CONFIRMED (nwdynam.ado's choice_coordination tertius()/four reproduce the direct-Mata-call crosscheck exactly through the real .ado command)"
