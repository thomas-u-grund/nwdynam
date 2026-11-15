* nwdynam plugin (native/dynam_sim.c) against the Mata reference evaluator of the choice sub-model (cscripts/ref_dynam_mata.do):
* random small event sequences with simultaneous events, event values (also zero and negative), the sequential option, windows
* on inertia and recip, weighted inertia/recip/indeg/outdeg, two-mode networks, opportunities, tie(), every subset of the 17 effects.
* Log likelihood and gradient must agree to 1e-9 (relative), the Hessian with the numerical derivative of the reference gradient to 1e-5.
local nwcore : env NWCORE
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do cscripts/ref_dynam_mata.do
do cscripts/test_nwdynam_plugin_lib.do
mata:
mata set matastrict off
void runchoice() {
	class DynamState scalar S
	class DynamNat scalar NT
	real scalar worst, worsth, c, n, N, ties, valued, seq, two, i, k, j, wi, wr, wti, wtr, wtd, wto, ll0, ll1, d, dh, nact, useopp, usetie
	real matrix ev, H1, hn, covmat, opp, tiem
	real rowvector active, th, g0, g1, tp, tm, gp, gm, windowvec, weightvec
	worst = 0; worsth = 0
	rseed(8812)
	for (c=1; c<=150; c++) {
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
		active = (runiform(1,17) :> (c <= 40 ? 0.8 : 0.5))
		if (sum(active)==0) active[3] = 1
		wi = (mod(c,4)==0 ? 1.5 : 1e300); wr = (mod(c,6)==0 ? 2 : 1e300)
		wti = (mod(c,8)==0); wtr = (mod(c,10)==0); wtd = (mod(c,12)==0); wto = (mod(c,13)==0)
		useopp = (mod(c,5)==2); usetie = (mod(c,4)==3 | active[17])
		active[17] = usetie
		covmat = (floor(runiform(1,n):*3) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n))
		opp = J(0,0,0)
		tiem = J(0,0,0)
		if (usetie) tiem = (runiform(n,n) :> 0.6)
		st_numscalar("__nwdynam_sequential", seq)
		S = DynamState(); S.init(ev, n)
		if (useopp) {
			opp = (runiform(S.nevents, n) :> 0.3)
			for (i=1;i<=S.nevents;i++) opp[i, S.events[i,2]] = 1
			if (two) opp = opp :* ((1..n) :> two)
			if (two) for (i=1;i<=S.nevents;i++) opp[i, S.events[i,2]] = 1
		}
		windowvec = (wi, wr)
		weightvec = (wti, wtr, wtd, wto)
		NT = DynamNat(); NT.open()
		NT.setup(1, &S, active, covmat, (wi, wr, 1e300, 1e300), weightvec, two, opp, tiem)
		nact = sum(active)
		th = (runiform(1,nact) :- 0.5) :* 0.5
		dynam_choice_loglik_grad_multi(th, active, covmat, windowvec, weightvec, two, opp, tiem, &S, ll0, g0)
		dynam_nat_eval(2, th, &NT, ll1, g1, H1)
		d = max(abs((ll0, g0) - (ll1, g1)) :/ (1 :+ abs((ll0, g0))))
		hn = J(cols(th), cols(th), 0)
		for (k=1;k<=cols(th);k++) {
			tp = th; tm = th; tp[k] = tp[k]+1e-5; tm[k] = tm[k]-1e-5
			dynam_choice_loglik_grad_multi(tp, active, covmat, windowvec, weightvec, two, opp, tiem, &S, j, gp)
			dynam_choice_loglik_grad_multi(tm, active, covmat, windowvec, weightvec, two, opp, tiem, &S, j, gm)
			hn[k,.] = (gp - gm) / 2e-5
		}
		dh = max(abs(hn - H1) :/ (1 :+ abs(hn)))
		worst = max((worst, d)); worsth = max((worsth, dh))
		NT.close()
		if (d > 1e-9 | dh > 1e-5) printf("case %g DIFF d=%g dh=%g active=%s ties=%g valued=%g seq=%g two=%g wi=%g wr=%g wt=%s opp=%g tie=%g\n", c, d, dh, invtokens(strofreal(active)), ties, valued, seq, two, wi, wr, invtokens(strofreal(weightvec)), useopp, usetie)
	}
	st_numscalar("worstll", worst); st_numscalar("worsth", worsth)
}
runchoice()
end
di "worst relative difference ll/grad: " scalar(worstll) "  Hessian vs numerical: " scalar(worsth)
assert scalar(worstll) < 1e-9
assert scalar(worsth) < 1e-5
di "PASS: plugin choice equals the Mata reference"
