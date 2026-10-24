* dev/dynam_unit20b_rateintercept_ado_crosscheck.do -- confirms
* nwdynam.ado's intercept option (batch 17, 2026-09-02) reproduces
* dev/dynam_unit20_rateintercept_crosscheck.do's own direct-Mata-call
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
import delimited "dev/dynam_crosscheck_intercept_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(intnet)

di as text "{hline 60}"
di as text "nwdynam intnet, submodel(rate) intercept indeg"
di as text "{hline 60}"
nwdynam intnet, submodel(rate) intercept indeg
matrix b0 = e(b)
scalar ll0 = e(ll)
di as text "e(effects) = `e(effects)'"
scalar d0a = abs(el(b0,1,1) - (-3.920622))
scalar d0b = abs(el(b0,1,2) - (-0.529166))
scalar d0c = abs(ll0 - (-51.8521))
di as text "abs diff: intercept=" %9.6f d0a "  indeg=" %9.6f d0b "  logLik=" %9.6f d0c
assert d0a < 0.001
assert d0b < 0.001
assert d0c < 0.001

di as text "{hline 60}"
di as text "Test: intercept rejected under submodel(choice)"
di as text "{hline 60}"
capture noisily nwdynam intnet, submodel(choice) intercept
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as text "{hline 60}"
di as text "Test: intercept rejected combined with indegwindow()"
di as text "{hline 60}"
capture noisily nwdynam intnet, submodel(rate) intercept indegwindow(3)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as text "{hline 60}"
di as text "Test: intercept + weightedindeg now fits successfully (batch 18,"
di as text "2026-09-02 - see dev/dynam_unit21_rateintercept_extras_crosscheck.do"
di as text "for the real cross-tool correctness check; this file's own scope"
di as text "stays batch 17's original reference numbers above)"
di as text "{hline 60}"
nwdynam intnet, submodel(rate) intercept weightedindeg
di as text "  fits successfully, e(effects)=`e(effects)'"
assert "`e(effects)'" == "Intercept indeg"

di as result "MATCH CONFIRMED (nwdynam.ado's intercept option reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
