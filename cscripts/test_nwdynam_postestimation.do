cscript

* Postestimation of nwdynam against goldfish 1.6.12 (R 4.6) and the R package changepoint 2.3: estat ic, estat eventlik, estat eventprob,
* estat outliers, estat changepoints. The files of dev/post (made by dev/dynam_post_crosscheck.R; the test does not need R) hold random event logs
* (choice, rate, rate with intercept, choice_coordination; composition change, time-varying covariates, windows, two-mode, opportunities,
* netterms) and goldfish's values: per-event log likelihoods and probabilities (1e-10), logLik, AIC, BIC, outlier indices and change points (exact).
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
capture adopath ++ `"`nwcore'/lib"'
adopath ++ `"`c(pwd)'"'
adopath ++ `"`c(pwd)'/lib"'
adopath ++ `"`c(pwd)'/lib/plugins/macos"'
global POSTONLY : env POSTONLY
global POSTNOPROB : env POSTNOPROB
global POSTDIR `"`c(pwd)'/dev/post"'
quietly do dev/dynam_post_crosscheck.do
prunall

* ---- 2. identities (choice: the event log likelihoods add up to the log likelihood, the probabilities of an event add up to 1), error paths
local nfail 0
capture noisily prun r06 choice "inertia same(g) sim(y) nodetrans" 1 ""
estat eventlik, frame(ev6)
assert abs(r(ll_events) - e(ll)) < 1e-8
assert r(n_cens) == 0
estat eventprob, frame(pr6)
frame pr6: collapse (sum) p, by(event)
frame pr6: assert abs(p - 1) < 1e-12
estat ic
assert r(N) == e(N)
assert abs(r(aic) - (-2 * e(ll) + 2 * colsof(e(b)))) < 1e-9
assert abs(r(bic) - (-2 * e(ll) + ln(e(N)) * colsof(e(b)))) < 1e-9
* estat outliers / changepoints: frames and returns
estat outliers, method(iqr) parameter(1) frame(ol6)
assert r(N) == e(N)
local nout = r(n_outliers)
frame ol6: quietly count if outlier
assert r(N) == `nout'
estat changepoints, window(3) frame(cp6)
assert r(window) == 3
assert colsof(r(cpts)) == r(n_cpts) + 1
assert el(r(cpts), 1, colsof(r(cpts))) == e(N)
* the table frame can be listed and carries goldfish's event numbers
frame ev6: assert event == _n
frame ev6: assert origevent >= 1
* error paths
capture estat eventlik, frame(ev6)
assert _rc == 110
estat eventlik, frame(ev6) replace
capture noisily estat eventprob
assert _rc == 198
capture estat eventprob, frame(pr6)
assert _rc == 110
capture estat eventprob, frame(pr6x) maxrows(10)
assert _rc != 0
capture estat foo
assert _rc == 321
quietly estat vce
assert rowsof(r(V)) == colsof(e(b))
capture estat outliers, method(xyz)
assert _rc == 198
capture estat outliers, window(-1)
assert _rc == 198
capture estat outliers, window(abc)
assert _rc == 198
capture estat outliers, method(hampel) window(1000)
assert _rc == 198
capture estat outliers, method(top) parameter(100000)
assert _rc == 198
capture estat changepoints, moment(skew)
assert _rc == 198
capture estat changepoints, method(xyz)
assert _rc == 198
capture estat changepoints, window(0)
assert _rc == 198
capture estat changepoints, window(2.5)
assert _rc == 198
capture estat changepoints, window(1000)
assert _rc == 198
capture estat changepoints, method(binseg) window(1000)
assert _rc == 198
* a later fit replaces the kept data; the estimates of the earlier fit are then stale
estimates store old
quietly prun r01 choice "inertia inertiawindow(4) recip recipwindow(14) indeg diff(y) present(present0) compchange(comp)" 1 ""
estimates restore old
capture estat ic
assert _rc == 498
quietly prun r01 choice "inertia inertiawindow(4) recip recipwindow(14) indeg diff(y) present(present0) compchange(comp)" 1 ""
estat ic
* the data of the fit is gone after clear all (the estimates stay)
frame change default
capture frame drop __nwdynam_fit
capture estat ic
assert _rc == 498
* no nwdynam estimates
ereturn clear
capture estat ic
assert _rc == 301
* a run of the commands leaves the data and the frames of the user as they were
quietly prun r01 choice "inertia inertiawindow(4) recip recipwindow(14) indeg diff(y) present(present0) compchange(comp)" 1 ""
quietly count
local n0 = r(N)
quietly estat ic
quietly estat outliers
quietly estat changepoints
quietly estat eventlik
quietly count
assert r(N) == `n0'
assert c(frame) == "default"
di as result "ALL TESTS PASSED (postestimation: identities and error paths)"
