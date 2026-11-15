* nwdynam plugin (native/dynam_sim.c) against the Mata reference evaluators (cscripts/ref_dynam_mata.do):
* rate sub-model without and with intercept, on random small event sequences with simultaneous events, event values (also
* zero and negative), the sequential option, windows, weighted degrees, two-mode networks, every subset of the five effects.
* Log likelihood and gradient must agree to 1e-9 (relative), the Hessian with the numerical derivative of the reference gradient
* to 1e-5, and for the model with intercept with the Hessian of the cached reference evaluator to 1e-9.
local nwcore : env NWCORE
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do cscripts/ref_dynam_mata.do
do cscripts/test_nwdynam_plugin_lib.do
mata:
mata set matastrict off
void runrate() {
	class DynamState scalar S
	class DynamNat scalar NT
	struct DynamRICache scalar C
	real scalar worst, worsth, worstc, c, n, N, ties, valued, seq, two, i, k, wi, wo, wgi, wgo, ll0, ll1, d, dh, dc, model, nact
	real matrix ev, H1, hn, covmat, H2
	real rowvector active, ego, tert, modvec, th, g0, g1, tp, tm, gp, gm, g2
	worst = 0; worsth = 0; worstc = 0
	rseed(4711)
	for (c=1; c<=120; c++) {
		model = 2 + mod(c,2)          // 2 = rate, 3 = rate with intercept
		n = 6 + runiformint(1,1,0,9); N = 30 + runiformint(1,1,0,80)
		ties = mod(c,3) != 1; valued = (mod(c,5)==0) + 2*(mod(c,11)==0); seq = (mod(c,7)==0)
		two = (mod(c,9)==0 ? 3 + runiformint(1,1,0,2) : 0)
		ev = mkev(n, N, ties, valued)
		if (two) {
			for (i=1;i<=N;i++) {
				ev[i,1] = runiformint(1,1,1,two)
				ev[i,2] = two + runiformint(1,1,1,n-two)
			}
		}
		active = (runiform(1,5) :> 0.4)
		if (two) {
			active[1] = 0
			active[3] = 0
			active[4] = 0
			active[5] = 0
		}
		if (sum(active)==0) active[2] = 1
		wi = (mod(c,4)==0 ? 1.5 : 1e300); wo = (mod(c,6)==0 ? 2 : 1e300)
		wgi = (mod(c,8)==0); wgo = (mod(c,10)==0)
		ego = runiform(1,n); tert = runiform(1,n)
		modvec = (wi, wo, wgi, wgo, two)
		st_numscalar("__nwdynam_sequential", seq); S = DynamState(); S.init(ev, n)
		covmat = (J(3,n,0) \ ego \ tert \ J(2,n,0))
		NT = DynamNat(); NT.open()
		NT.setup(model, &S, active, covmat, (1e300,1e300,wi,wo), (0,0,wgi,wgo), two, J(0,0,0), J(0,0,0))
		nact = sum(active)
		if (model == 3) {
			th = (runiform(1,1+nact) :- 0.5) :* 0.6
			dynam_rateint_loglik_grad_multi(th, active, ego, tert, modvec, &S, ll0, g0)
			dynam_nat_eval(2, th, &NT, ll1, g1, H1)
			C = dynam_rateint_build_cache(&S, active, ego, tert, modvec)
			dynam_rateint_cached_eval(2, th, &C, 0, d, g2, H2)
			dc = max(abs(H1 - H2) :/ (1 :+ abs(H2)))
			worstc = max((worstc, dc))
		}
		else {
			th = (runiform(1,nact) :- 0.5) :* 1.2
			dynam_rate_loglik_grad_multi(th, active, ego, tert, wi, wo, wgi, wgo, two, &S, ll0, g0)
			dynam_nat_eval(2, th, &NT, ll1, g1, H1)
			dc = 0
		}
		d = max(abs((ll0, g0) - (ll1, g1)) :/ (1 :+ abs((ll0, g0))))
		// numerical Hessian from the reference gradient
		hn = J(cols(th), cols(th), 0)
		for (k=1;k<=cols(th);k++) {
			tp = th; tm = th; tp[k] = tp[k]+1e-5; tm[k] = tm[k]-1e-5
			if (model == 3) {
				dynam_rateint_loglik_grad_multi(tp, active, ego, tert, modvec, &S, i, gp)
				dynam_rateint_loglik_grad_multi(tm, active, ego, tert, modvec, &S, i, gm)
			}
			else {
				dynam_rate_loglik_grad_multi(tp, active, ego, tert, wi, wo, wgi, wgo, two, &S, i, gp)
				dynam_rate_loglik_grad_multi(tm, active, ego, tert, wi, wo, wgi, wgo, two, &S, i, gm)
			}
			hn[k,.] = (gp - gm) / 2e-5
		}
		dh = max(abs(hn - H1) :/ (1 :+ abs(hn)))
		worst = max((worst, d)); worsth = max((worsth, dh))
		NT.close()
		if (d > 1e-9 | dh > 1e-5 | dc > 1e-9) printf("case %g DIFF model=%g d=%g dh=%g dc=%g active=%s ties=%g valued=%g seq=%g two=%g wi=%g wo=%g wgi=%g wgo=%g\n", c, model, d, dh, dc, invtokens(strofreal(active)), ties, valued, seq, two, wi, wo, wgi, wgo)
	}
	st_numscalar("worstll", worst); st_numscalar("worsth", worsth); st_numscalar("worstc", worstc)
}
runrate()
end
di "worst relative difference ll/grad: " scalar(worstll) "  Hessian vs numerical: " scalar(worsth) "  vs cached: " scalar(worstc)
assert scalar(worstll) < 1e-9
assert scalar(worsth) < 1e-5
assert scalar(worstc) < 1e-9
di "PASS: plugin rate / rate with intercept equal the Mata reference"
