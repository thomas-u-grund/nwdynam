* nwdynam rate sub-model with intercept: speed and converged estimates.
* The fit on the manufacturing emails (first 14 days, 2,738 events, 129 actors) finishes in seconds and has the converged estimates
* (the plugin against the Mata reference evaluators on random cases: test_nwdynam_plugin_rate.do).
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do nwdynam.ado

* manufacturing emails, first 14 days: seconds instead of an hour, converged estimates
use `"`nwcore'/data/manufacturing.dta"', clear
keep if t < 14
assert _N == 2738
nwset sender receiver, eventtime(t) name(mail)
timer clear 77
timer on 77
nwdynam mail, submodel(rate) intercept indeg
timer off 77
quietly timer list 77
di as text "fit time " r(t77) " s"
assert r(t77) < 120
matrix b = e(b)
assert abs(b[1,1] - 0.573097769813) < 1e-7
assert abs(b[1,2] - 0.021186290096) < 1e-7
assert abs(e(ll) - (-703.7556798644)) < 1e-6
assert e(nodes) == 129
di as text "PASS intercept rate fit on 2,738 manufacturing events"
di as text "ALL PASS: test_nwdynam_rateint_fast"
