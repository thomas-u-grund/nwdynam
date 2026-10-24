* dev/dynam_unit13b_opportunities_ado_crosscheck.do -- confirms the
* nwdynam.ado opportunities() option (batch 10, 2026-09-02) reproduces
* dev/dynam_unit13_opportunities_crosscheck.do's own direct-Mata-call
* numbers when run through the REAL .ado command, on the SAME data
* (dev/dynam_crosscheck_oppl_events.csv / _availability.csv, built by
* dev/dynam_unit13_opportunities_crosscheck.R against real goldfish
* 1.6.12 - see that file's own header comment for the reference
* numbers and mechanics discovery).

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
import delimited "dev/dynam_crosscheck_oppl_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(oppnet)

di as text "{hline 60}"
di as text "nwdynam oppnet, submodel(choice) inertia indeg  (no opportunities())"
di as text "{hline 60}"
nwdynam oppnet, submodel(choice) inertia indeg
matrix b0 = e(b)
scalar ll0 = e(ll)
scalar d0a = abs(el(b0,1,1) - 0.45766)
scalar d0b = abs(el(b0,1,2) - (-0.54225))
scalar d0c = abs(ll0 - (-56.75))
di as text "abs diff: inertia=" %9.6f d0a "  indeg=" %9.6f d0b "  logLik=" %9.4f d0c
assert d0a < 0.001
assert d0b < 0.001
assert d0c < 0.01

di as text "{hline 60}"
di as text "nwdynam oppnet, submodel(choice) inertia indeg opportunities(event actor)"
di as text "{hline 60}"
import delimited "dev/dynam_crosscheck_oppl_availability.csv", clear varnames(1)
rename event oppevent
rename actor oppactor
nwdynam oppnet, submodel(choice) inertia indeg opportunities(oppevent oppactor)
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - 0.52986)
scalar d1b = abs(el(b1,1,2) - (-0.52236))
scalar d1c = abs(ll1 - (-52.104))
di as text "abs diff: inertia=" %9.6f d1a "  indeg=" %9.6f d1b "  logLik=" %9.4f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.01

di as text "{hline 60}"
di as text "opportunities() rejected under submodel(rate)"
di as text "{hline 60}"
capture noisily nwdynam oppnet, submodel(rate) opportunities(oppevent oppactor)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as result "MATCH CONFIRMED (nwdynam.ado's opportunities() reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
