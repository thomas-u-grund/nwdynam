* dev/dynam_unit21b_rateintercept_extras_ado_crosscheck.do -- confirms
* nwdynam.ado's intercept option, now combined with two-mode and
* weighted=TRUE (batch 18), reproduces dev/dynam_unit21_
* rateintercept_extras_crosscheck.do's own direct-Mata-call numbers
* through the REAL .ado command, on the SAME data. Also confirms
* indegwindow()/outdegwindow() combined with intercept is still rejected
* (a real, disclosed architectural gap - see nwdynam.ado's own intercept
* validation block).

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

di as text "{hline 60}"
di as text "Sub-check 1: intercept + outdeg, two-mode, through nwdynam"
di as text "{hline 60}"
nwclear
import delimited "dev/dynam_crosscheck_twomode_intercept_events.csv", clear varnames(1)
encode sender, gen(sendernum)
encode receiver, gen(receivernum)
nwset sendernum receivernum, twomode eventtime(time) name(twomodeintnet)
nwdynam twomodeintnet, submodel(rate) intercept outdeg
di as text "e(effects) = `e(effects)'"
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - (-3.52957420))
scalar d1b = abs(el(b1,1,2) - (0.00931169))
scalar d1c = abs(ll1 - (-54.25246354))
di as text "abs diff: intercept=" %9.6f d1a "  outdeg=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001
di as result "  SUB-CHECK 1 PASSED"

di as text "{hline 60}"
di as text "Sub-check 2: intercept + outdeg WEIGHTED, through nwdynam"
di as text "{hline 60}"
nwclear
import delimited "dev/dynam_crosscheck_rateintweighted_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(weightedintnet)
nwdynam weightedintnet, submodel(rate) intercept weightedoutdeg
di as text "e(effects) = `e(effects)'"
matrix b2 = e(b)
scalar ll2 = e(ll)
scalar d2a = abs(el(b2,1,1) - (-2.90509662))
scalar d2b = abs(el(b2,1,2) - (-0.63982127))
scalar d2c = abs(ll2 - (-53.25960504))
di as text "abs diff: intercept=" %9.6f d2a "  outdegweighted=" %9.6f d2b "  logLik=" %9.6f d2c
assert d2a < 0.001
assert d2b < 0.001
assert d2c < 0.001
di as result "  SUB-CHECK 2 PASSED"

di as text "{hline 60}"
di as text "Sub-check 3: intercept + outdegwindow() is still rejected"
di as text "{hline 60}"
capture noisily nwdynam weightedintnet, submodel(rate) intercept outdegwindow(10)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as result "MATCH CONFIRMED (nwdynam.ado's intercept option reproduces the direct-Mata-call crosscheck exactly, combined with two-mode and weighted=TRUE, through the real .ado command; window() combination correctly still rejected)"
