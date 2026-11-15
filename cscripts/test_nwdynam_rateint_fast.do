* nwdynam rate sub-model with intercept: the fast evaluators (cached Mata, native kernel) against the plain Mata evaluator.
* 1. Over 60 random event sequences (simultaneous events, event values, sequential option, windows, weighted degrees,
*    two-mode, every subset of the five effects) the cached evaluator reproduces log likelihood and gradient of the plain
*    evaluator, and its Hessian matches the numerical derivative of that gradient.
* 2. The fit on the manufacturing emails (first 14 days, 2,738 events) finishes in seconds and has the converged estimates.
* 3. With the native plugin present, the native fit and the cached fit agree.
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do nwdynam.ado
mata:
mata set matastrict off
real matrix mkev(real scalar n, real scalar N, real scalar ties, real scalar valued) {
	real matrix ev
	real scalar i, a, b, t, v
	ev = J(N,4,1); t = 0
	for (i=1;i<=N;i++) {
		a = runiformint(1,1,1,n); b = runiformint(1,1,1,n); while (b==a) b = runiformint(1,1,1,n)
		if (!(ties & runiform(1,1) < 0.3)) t = t + 0.1 + runiform(1,1)*2
		v = 1
		if (valued) v = 1 + runiformint(1,1,0,3)
		ev[i,.] = (a,b,t,v)
	}
	return(ev)
}
void runcmp() {
	class DynamState scalar S
	struct DynamRICache scalar C
	real scalar worst, worsth, c, n, N, ties, valued, seq, two, i, k, wi, wo, wgi, wgo, ll0, ll1, d, dh, a1, a2
	real matrix ev, H1, hn
	real rowvector active, ego, tert, modvec, th, g0, g1, tp, tm, gp, gm
worst = 0
worsth = 0
rseed(12345)
for (c=1; c<=60; c++) {
	n = 8 + runiformint(1,1,0,8); N = 40 + runiformint(1,1,0,60)
	ties = mod(c,2); valued = (mod(c,3)==0); seq = (mod(c,5)==0)
	two = (mod(c,7)==0 ? 3 + runiformint(1,1,0,2) : 0)
	ev = mkev(n, N, ties, valued)
	if (two) {
		for (i=1;i<=N;i++) {
			ev[i,1] = runiformint(1,1,1,two)
			ev[i,2] = two + runiformint(1,1,1,n-two)
		}
	}
	active = (runiform(1,5) :> 0.5)
	if (two) active[1] = 0
	if (two) active[3] = 0
	if (two) active[4] = 0
	if (two) active[5] = 0
	if (sum(active)==0) active[2] = 1
	wi = (mod(c,4)==0 ? 1.5 : 1e300); wo = (mod(c,6)==0 ? 2 : 1e300)
	wgi = (mod(c,9)==0); wgo = (mod(c,8)==0)
	ego = runiform(1,n); tert = runiform(1,n)
	modvec = (wi, wo, wgi, wgo, two)
	S = DynamState(); S.sequential = seq; S.init(ev, n)
	th = (runiform(1,1+sum(active)) :- 0.5) :* 0.6
	dynam_rateint_loglik_grad_multi(th, active, ego, tert, modvec, &S, ll0, g0)
	C = dynam_rateint_build_cache(&S, active, ego, tert, modvec)
	dynam_rateint_cached_eval(2, th, &C, 0, ll1, g1, H1)
	d = max(abs((ll0, g0) - (ll1, g1)) :/ (1 :+ abs((ll0, g0))))
	// numerical check of the Hessian by differencing the old gradient
	hn = J(cols(th), cols(th), 0)
	for (k=1;k<=cols(th);k++) {
		tp = th; tm = th; tp[k] = tp[k]+1e-5; tm[k] = tm[k]-1e-5
		dynam_rateint_loglik_grad_multi(tp, active, ego, tert, modvec, &S, a1, gp)
		dynam_rateint_loglik_grad_multi(tm, active, ego, tert, modvec, &S, a2, gm)
		hn[k,.] = (gp - gm) / 2e-5
	}
	dh = max(abs(hn - H1) :/ (1 :+ abs(hn)))
	worst = max((worst, d)); worsth = max((worsth, dh))
	if (d > 1e-9 | dh > 1e-5) printf("case %g DIFF d=%g dh=%g active=%s ties=%g valued=%g seq=%g two=%g\n", c, d, dh, invtokens(strofreal(active)), ties, valued, seq, two)
}
st_numscalar("worstll", worst)
	st_numscalar("worsth", worsth)
}
runcmp()

void natcmp() {
	class DynamState scalar S
	real matrix ev
	real scalar c, ok, nat
	real rowvector modvec, act
	nat = 0
	st_numscalar("natused", 0)
	st_numscalar("natworst", 0)
	for (c=1; c<=6; c++) {
		ev = mkev(10 + c, 120 + 20*c, mod(c,2), (c > 4))
		S = DynamState()
		S.init(ev, 10 + c)
		modvec = (1e300, 1e300, 0, 0, 0)
		act = (1, mod(c,3) != 0, 0, 0, 0)
		if (dynam_rateint_native_fit(&S, act, modvec, "bn", "Vn", "lln")) {
			nat = nat + 1
			ok = dynam_rateint_cached_fit(&S, act, J(1,10+c,0), J(1,10+c,0), modvec, "bc", "Vc", "llc")
			st_numscalar("natworst", max((st_numscalar("natworst"), max(abs(st_matrix("bn") - st_matrix("bc")) :/ (1 :+ abs(st_matrix("bc")))), max(abs(st_matrix("Vn") - st_matrix("Vc")) :/ (1 :+ abs(st_matrix("Vc")))), abs(st_numscalar("lln") - st_numscalar("llc")) / (1 + abs(st_numscalar("llc"))), 1 - ok)))
		}
	}
	st_numscalar("natused", nat)
}
natcmp()
end

di as text "worst relative difference (log likelihood, gradient): " scalar(worstll) ", Hessian against numerical derivative: " scalar(worsth)
assert scalar(worstll) < 1e-9
assert scalar(worsth) < 1e-5
di as text "PASS cached evaluator equals the plain evaluator on 60 random cases"

* native kernel against the cached fit (only where a plugin for this platform exists)
if scalar(natused) > 0 {
	di as text "native fits compared: " scalar(natused) ", worst relative difference " scalar(natworst)
	assert scalar(natworst) < 1e-7
	di as text "PASS native kernel equals the cached fit"
}
else di as text "no native plugin on this platform: native comparison skipped"

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
