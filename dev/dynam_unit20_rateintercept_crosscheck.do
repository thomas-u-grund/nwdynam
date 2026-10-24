* dev/dynam_unit20_rateintercept_crosscheck.do -- WITH-INTERCEPT rate
* sub-model, batch 17, resolved 2026-09-02 - the genuinely
* continuous-time competing-risks hazard variant, against real
* goldfish subset fits on a 6-actor/10-event toy directed network with
* real (unevenly-spaced) timestamps - direct Mata-call discipline.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit20_rateintercept_crosscheck.R for the R side):
*   rate ~ 1 + indeg:   intercept=-3.920622, indeg=-0.529166, logLik=-51.8521
*   rate ~ 1 + outdeg:  intercept=-4.4196085, outdeg=0.2009892, logLik=-52.18619

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
do unw_dynam.do

import delimited "dev/dynam_crosscheck_intercept_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat = (sender, receiver, eventtime)
n_actors = 6
zero6 = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "rate ~ 1 + indeg"
di as text "{hline 60}"
mata: DynamRateInterceptFitMulti(eventmat, n_actors, (1,0,0,0,0), zero6, zero6, 1e300, 1e300, 0, 0, 0, "b0", "V0", "ll0")
assert colsof(b0) == 2
di as text "intercept=" %9.6f el(b0,1,1) "  indeg=" %9.6f el(b0,1,2) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - (-3.920622))
scalar d0b = abs(el(b0,1,2) - (-0.529166))
scalar d0c = abs(ll0 - (-51.8521))
di as text "abs diff: intercept=" %9.6f d0a "  indeg=" %9.6f d0b "  logLik=" %9.6f d0c
assert d0a < 0.001
assert d0b < 0.001
assert d0c < 0.001

di as text "{hline 60}"
di as text "rate ~ 1 + outdeg"
di as text "{hline 60}"
mata: DynamRateInterceptFitMulti(eventmat, n_actors, (0,1,0,0,0), zero6, zero6, 1e300, 1e300, 0, 0, 0, "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "intercept=" %9.6f el(b1,1,1) "  outdeg=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - (-4.4196085))
scalar d1b = abs(el(b1,1,2) - 0.2009892)
scalar d1c = abs(ll1 - (-52.18619))
di as text "abs diff: intercept=" %9.6f d1a "  outdeg=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as result "MATCH CONFIRMED (WITH-INTERCEPT rate sub-model reproduces real goldfish exactly on both indeg and outdeg)"
