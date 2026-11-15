* nwdynam plugin (native/dynam_sim.c) against the Mata reference evaluator of choice_coordination (cscripts/ref_dynam_mata.do):
* random small undirected event sequences with simultaneous events, event values (also zero and negative), the sequential option,
* every subset of the 12 effects, tie(). Log likelihood and gradient must agree to 1e-9 (relative).
local nwcore : env NWCORE
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do cscripts/ref_dynam_mata.do
do cscripts/test_nwdynam_plugin_lib.do
mata:
mata set matastrict off
void runcoord() {
	class DynamState scalar S
	class DynamNat scalar NT
	real scalar worst, c, n, N, ties, valued, seq, d, nact, llm, lln
	real matrix ev, covmat, tiem, H
	real rowvector active, th, gm, gn
	worst = 0
	rseed(31337)
	for (c=1; c<=100; c++) {
		n = 5 + runiformint(1,1,0,6); N = 25 + runiformint(1,1,0,40)
		ties = mod(c,3) != 1; valued = (mod(c,5)==0) + 2*(mod(c,11)==0); seq = (mod(c,7)==0)
		ev = mkev(n, N, ties, valued)
		active = (runiform(1,12) :> (c <= 30 ? 0.8 : 0.5))
		if (sum(active)==0) active[1] = 1
		covmat = (floor(runiform(1,n):*3) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n) \ runiform(1,n))
		tiem = J(0,0,0)
		if (active[12]) tiem = (runiform(n,n) :> 0.6)
		st_numscalar("__nwdynam_sequential", seq)
		S = DynamState(); S.init(ev, n)
		NT = DynamNat(); NT.open()
		NT.setup(4, &S, active, covmat, (1e300,1e300,1e300,1e300), (0,0,0,0), 0, J(0,0,0), tiem)
		nact = sum(active)
		th = (runiform(1,nact) :- 0.5) :* 0.8
		dynam_coord_loglik_grad_multi(th, active, covmat, tiem, &S, llm, gm)
		dynam_nat_eval(1, th, &NT, lln, gn, H)
		d = max(abs((llm, gm) - (lln, gn)) :/ (1 :+ abs((llm, gm))))
		worst = max((worst, d))
		NT.close()
		if (d > 1e-9) printf("case %g DIFF d=%g active=%s ties=%g valued=%g seq=%g n=%g\n", c, d, invtokens(strofreal(active)), ties, valued, seq, n)
	}
	st_numscalar("worstll", worst)
}
runcoord()
end
di "worst relative difference ll/grad: " scalar(worstll)
assert scalar(worstll) < 1e-9
di "PASS: plugin choice_coordination equals the Mata reference"
