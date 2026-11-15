cscript

* nwdynam plugin on a small hand-built event stream and on a random stream, against the Mata reference evaluators of the unit-1
* models (cscripts/ref_dynam_mata.do):
*  (1) the plugin is found and answers with the protocol of this package;
*  (2) log likelihood and gradient of the three-effect choice model and the two-effect rate model, several theta vectors;
*  (3) full fits of nwdynam: the estimates equal the maximum found by Mata's optimize() on the reference evaluator.
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do cscripts/ref_dynam_mata.do

mata: mata set matastrict off
mata: printf("DynamNativeAvailable() = %g (path: %s)\n", DynamNativeAvailable(), DynamNativePluginPath())
mata: assert(DynamNativeAvailable() == 1)

mata:
void cmpunit1() {
	class DynamState scalar S
	class DynamNat scalar NT
	real matrix thetas_c, thetas_r, H
	real rowvector theta, gm, gn
	real scalar row, llm, lln
	S = DynamState()
	S.init((1,2,1 \ 2,1,2 \ 3,1,3 \ 1,2,4 \ 3,1,5 \ 4,3,6 \ 2,4,7 \ 5,1,8), 5)
	thetas_c = (0,0,0 \ 0.5,-0.3,0.2 \ -1.1,2.4,-0.6 \ 3,3,3)
	NT = DynamNat()
	NT.open()
	NT.setup(1, &S, (1,1,1,J(1,14,0)), J(7,5,0), (1e300,1e300,1e300,1e300), (0,0,0,0), 0, J(0,0,0), J(0,0,0))
	for (row=1; row<=rows(thetas_c); row++) {
		theta = thetas_c[row,.]
		dynam_choice_loglik_grad_unit1(theta, &S, llm, gm)
		dynam_nat_eval(1, theta, &NT, lln, gn, H)
		printf("choice row %g: |ll diff|=%g  |grad diff|=%g\n", row, abs(llm - lln), max(abs(gm - gn)))
		assert(abs(llm - lln) < 1e-9)
		assert(max(abs(gm - gn)) < 1e-9)
	}
	NT.close()
	thetas_r = (0,0 \ 0.5,-0.3 \ -1.1,2.4 \ 3,3)
	NT = DynamNat()
	NT.open()
	NT.setup(2, &S, (1,1,0,0,0), J(7,5,0), (1e300,1e300,1e300,1e300), (0,0,0,0), 0, J(0,0,0), J(0,0,0))
	for (row=1; row<=rows(thetas_r); row++) {
		theta = thetas_r[row,.]
		dynam_rate_loglik_grad_unit1(theta, &S, llm, gm)
		dynam_nat_eval(1, theta, &NT, lln, gn, H)
		printf("rate row %g: |ll diff|=%g  |grad diff|=%g\n", row, abs(llm - lln), max(abs(gm - gn)))
		assert(abs(llm - lln) < 1e-9)
		assert(max(abs(gm - gn)) < 1e-9)
	}
	NT.close()
}
cmpunit1()

// maximum of the reference likelihood with Mata's optimize() (BFGS, numerical Hessian)
void refmax(real matrix ev, real scalar n, real scalar choice, string scalar bname, string scalar llname) {
	class DynamState scalar S
	transmorphic S_opt
	S = DynamState()
	S.init(ev, n)
	S_opt = optimize_init()
	optimize_init_evaluator(S_opt, choice ? &dynam_choice_eval_unit1() : &dynam_rate_eval_unit1())
	optimize_init_evaluatortype(S_opt, "d1")
	optimize_init_params(S_opt, J(1, choice ? 3 : 2, 0.01))
	optimize_init_argument(S_opt, 1, &S)
	optimize_init_technique(S_opt, "bfgs")
	optimize_init_which(S_opt, "max")
	optimize_init_tracelevel(S_opt, "none")
	optimize_init_conv_warning(S_opt, "off")
	optimize_init_conv_ptol(S_opt, 1e-12)
	optimize_init_conv_vtol(S_opt, 1e-14)
	optimize_init_conv_nrtol(S_opt, 1e-12)
	optimize_init_conv_maxiter(S_opt, 200)
	(void) optimize(S_opt)
	st_matrix(bname, optimize_result_params(S_opt))
	st_numscalar(llname, optimize_result_value(S_opt))
}
end

nwclear
clear
set obs 200
gen long sender = .
gen long receiver = .
gen eventtime = _n
mata:
rseed(20260902)
n = 8
for (i=1; i<=200; i++) {
	pair = runiformint(1, 2, 1, n)
	while (pair[1] == pair[2]) pair = runiformint(1, 2, 1, n)
	st_store(i, "sender", pair[1])
	st_store(i, "receiver", pair[2])
}
end
nwset sender receiver, eventtime(eventtime) name(cmpnet)
_nwsyntax cmpnet
mata: __ev = *(`netobj'->get_eventlist())

nwdynam cmpnet
matrix b_nat = e(b)
scalar ll_nat = e(ll)
mata: refmax(__ev, 8, 1, "b_ref", "ll_ref")
mata: assert(max(abs(st_matrix("b_nat") - st_matrix("b_ref"))) < 1e-5)
mata: assert(abs(st_numscalar("ll_nat") - st_numscalar("ll_ref")) < 1e-8)
di as text "choice: plugin fit " b_nat[1,1] " " b_nat[1,2] " " b_nat[1,3] "   reference " b_ref[1,1] " " b_ref[1,2] " " b_ref[1,3]

nwdynam cmpnet, submodel(rate)
matrix b_nat = e(b)
scalar ll_nat = e(ll)
mata: refmax(__ev, 8, 0, "b_ref", "ll_ref")
mata: assert(max(abs(st_matrix("b_nat") - st_matrix("b_ref"))) < 1e-5)
mata: assert(abs(st_numscalar("ll_nat") - st_numscalar("ll_ref")) < 1e-8)
di as text "rate:   plugin fit " b_nat[1,1] " " b_nat[1,2] "   reference " b_ref[1,1] " " b_ref[1,2]
mata: mata drop __ev

di "test_nwdynam_native: ALL OK"
