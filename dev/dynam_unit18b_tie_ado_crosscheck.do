* dev/dynam_unit18b_tie_ado_crosscheck.do -- confirms nwdynam.ado's
* tie() option (cross-network effects, v1 scope, 2026-09-02) reproduces
* dev/dynam_unit18_tie_crosscheck.do's own direct-Mata-call numbers
* through the REAL .ado command, on the SAME data.

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
import delimited "dev/dynam_crosscheck_tie_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
nwset senderidx receiveridx, eventtime(time) name(tienet)

* Load the SAME static matrix as an ordinary (non-event) network via
* nwset's own mat() option - the package's own established
* "declare a network from a Stata matrix" convention.
preserve
import delimited "dev/dynam_crosscheck_tie_staticmat.csv", clear varnames(1)
mkmat v1-v6, matrix(staticmat)
restore
nwset, mat(staticmat) name(staticnet) directed

di as text "{hline 60}"
di as text "nwdynam tienet, submodel(choice) tie(staticnet)"
di as text "{hline 60}"
nwdynam tienet, submodel(choice) tie(staticnet)
matrix b0 = e(b)
scalar ll0 = e(ll)
scalar d0a = abs(el(b0,1,1) - 1.791714)
scalar d0b = abs(ll0 - (-17.99618))
di as text "abs diff: tie=" %9.6f d0a "  logLik=" %9.6f d0b
assert d0a < 0.001
assert d0b < 0.001

di as text "{hline 60}"
di as text "nwdynam tienet, submodel(choice) inertia tie(staticnet)"
di as text "{hline 60}"
nwdynam tienet, submodel(choice) inertia tie(staticnet)
matrix b1 = e(b)
scalar ll1 = e(ll)
scalar d1a = abs(el(b1,1,1) - (-2.789568))
scalar d1b = abs(el(b1,1,2) - 2.287039)
scalar d1c = abs(ll1 - (-15.40548))
di as text "abs diff: inertia=" %9.6f d1a "  tie=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as text "{hline 60}"
di as text "Test: tie() rejected under submodel(rate)"
di as text "{hline 60}"
capture noisily nwdynam tienet, submodel(rate) tie(staticnet)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as text "{hline 60}"
di as text "Test: tie() with mismatched actor count rejected"
di as text "{hline 60}"
nwrandom 5, prob(.3) name(smallnet)
capture noisily nwdynam tienet, submodel(choice) tie(smallnet)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as text "{hline 60}"
di as text "Test: tie() with an event-type (temporal) network rejected"
di as text "{hline 60}"
capture noisily nwdynam tienet, submodel(choice) tie(tienet)
assert _rc == 198
di as text "  correctly rejected (rc=" _rc ")"

di as result "MATCH CONFIRMED (nwdynam.ado's tie() reproduces the direct-Mata-call crosscheck exactly through the real .ado command)"
