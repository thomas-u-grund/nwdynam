* dev/dynam_unit19b_coordination_tie_ado_crosscheck.do -- confirms
* nwdynam.ado's submodel(choice_coordination) tie() option (batch 16,
* 2026-09-02) reproduces dev/dynam_unit19_coordination_tie_
* crosscheck.do's own direct-Mata-call numbers through the REAL .ado
* command, on the SAME data.

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
import delimited "dev/dynam_crosscheck_coordtie_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(coordtienet) undirected

preserve
import delimited "dev/dynam_crosscheck_coordtie_staticmat.csv", clear varnames(1)
mkmat v1-v6, matrix(staticmat)
restore
nwset, mat(staticmat) name(statnet6) directed

di as text "{hline 60}"
di as text "nwdynam coordtienet, submodel(choice_coordination) inertia tie(statnet6)"
di as text "{hline 60}"
nwdynam coordtienet, submodel(choice_coordination) inertia tie(statnet6)
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - (-0.3867456))
scalar d1b = abs(el(b1,1,2) - 0.6756129)
scalar d1c = abs(ll1 - (-31.43784))
di as text "abs diff: inertia=" %9.6f d1a "  tie=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as result "MATCH CONFIRMED (nwdynam.ado's choice_coordination tie() reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
