* dev/dynam_unit15b_coordination_closure_ado_crosscheck.do -- confirms
* nwdynam.ado's submodel(choice_coordination) nodetrans/trans options
* (batch 12, 2026-09-02) reproduce dev/dynam_unit15_coordination_
* closure_crosscheck.do's own direct-Mata-call numbers through the
* REAL .ado command, on the SAME data.

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

di as text "{hline 60}"
di as text "nwdynam coordnet3, submodel(choice_coordination) inertia nodetrans trans"
di as text "{hline 60}"
nwdynam coordnet3, submodel(choice_coordination) inertia nodetrans trans
matrix b2 = e(b)
scalar ll2 = e(ll)
scalar d2a = abs(el(b2,1,1) - (-0.7852777))
scalar d2b = abs(el(b2,1,2) - 0.2335562)
scalar d2c = abs(el(b2,1,3) - (-0.1884200))
scalar d2d = abs(ll2 - (-31.5738))
di as text "abs diff: inertia=" %9.6f d2a "  nodeTrans=" %9.6f d2b "  trans=" %9.6f d2c "  logLik=" %9.6f d2d
assert d2a < 0.001
assert d2b < 0.001
assert d2c < 0.001
assert d2d < 0.001

di as result "MATCH CONFIRMED (nwdynam.ado's choice_coordination nodetrans/trans reproduce the direct-Mata-call crosscheck exactly through the real .ado command)"
