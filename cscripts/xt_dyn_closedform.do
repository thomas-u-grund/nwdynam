cscript

* xt_dyn_closedform.do -- DRAFT, not yet run.
* nwdynam against numbers derived by hand: scores and log likelihood at theta = 0 on a 4-event sequence (sequential
* and simultaneous), and closed-form maximum likelihood estimates where an effect splits every risk set into two
* groups. All values were also checked by brute force / numerical maximisation in R (xt_NOTES.md).
*
* Closed forms (W40 toy sequence, 4 actors, 40 events; x01 = (0,0,1,1), x4 = (0,0,0,1), next = (2,3,4,1)):
*   rate ego(x01):  a = 17 events sent by actors 3, 4; b = 23; n1 = n0 = 2:
*                   theta = ln(a n0/(b n1)), ll = a ln(a/(N n1)) + b ln(b/(N n0)), se = sqrt(N/(a b))
*   rate intercept ego(x01): two Poisson rates over T = t_40 - t_1 = 39:
*                   Intercept = ln(b/(n0 T)), ego = ln(a n0/(b n1)), ll = b ln(b/(n0 T)) + a ln(a/(n1 T)) - N,
*                   se(Intercept) = sqrt(1/b), se(ego) = sqrt(1/a + 1/b)
*   choice alter(x4): the 30 events not sent by actor 4 have one candidate with x4 = 1 among 3; a = 7 of them choose
*                   actor 4, b = 23 do not; the 10 events of actor 4 have a constant risk set (ln(1/3) each):
*                   theta = ln(2a/b), ll = a ln(a/30) + b ln(b/60) + 10 ln(1/3), se = sqrt(30/(a b))
*   choice tie(next): one candidate with tie = 1 among 3 at every event; a = 6 events go to next(sender):
*                   theta = ln(2a/b), ll = a ln(a/40) + b ln(b/80), se = sqrt(40/(a b))
* Expected runtime: about 20 seconds (Stata 17 BE).
* Run from the repository root; the test runner sets the adopath. No sysdir/adopath here.

capture program drop nwinstall
local nwcore : env NWCORE
if `"`nwcore'"' != "" quietly do `"`nwcore'/unw_core.do"'
quietly do unw_dynam.do

mata: XT_W40 = (1,3 \ 1,3 \ 1,3 \ 1,3 \ 2,1 \ 4,2 \ 2,4 \ 2,4 \ 1,3 \ 4,3 \ 3,1 \ 4,2 \ 3,2 \ 4,2 \ 1,3 \ 1,3 \ 2,4 \ 4,2 \ 2,4 \ 2,4 \ 3,1 \ 3,1 \ 1,2 \ 1,3 \ 3,4 \ 4,1 \ 1,3 \ 3,1 \ 1,2 \ 1,3 \ 4,2 \ 2,1 \ 4,2 \ 4,2 \ 3,2 \ 2,3 \ 1,3 \ 2,4 \ 4,2 \ 2,3)

capture program drop xt_w40
program xt_w40
	clear
	quietly set obs 40
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store(., ("sender", "receiver"), XT_W40)
	quietly gen double t = _n
end

capture program drop xt_cf
program xt_cf
	* xt_cf "label" ll b1 se1 [b2 se2 ...]: closed-form values; b to 1e-4, se to 1e-3 relative, ll to 1e-8 relative
	gettoken label 0 : 0
	gettoken ll 0 : 0
	local k = colsof(e(b))
	forvalues j = 1/`k' {
		gettoken b 0 : 0
		gettoken se 0 : 0
		assert abs(el(e(b), 1, `j') - (`b')) < 1e-4 * (1 + abs(`b'))
		assert reldif(sqrt(el(e(V), `j', `j')), `se') < 1e-3
	}
	assert reldif(e(ll), `ll') < 1e-8
	assert "`e(unidentified)'" == ""
	di as text "  PASS: `label' (b " %10.6f el(e(b),1,1) ", ll " %12.6f e(ll) ")"
end

capture program drop xt_cov
program xt_cov
	args net
	nwload `net', xvars
	capture assert _nwnode == "n" + string(_n)
	if _rc di as text "  note: actor labels are not n1..n" _N " - check the actor order used below"
end

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: scores at theta = 0, events 1->2, 2->1, 1->2, 1->3 on 3 actors (hand-derived)"
di as text "{hline 60}"
* times (1,2,3,4), sequential: choice (inertia, recip, indeg) = (0, 1/2, 0), rate (indeg, outdeg) = (4/3, 1/3)
* times (1,2,2,3), simultaneous: events 2 and 3 both see the state after event 1:
*                   choice = (0, 0, 0), rate = (2/3, 2/3)
* ll(0): choice -4 ln 2 (2 candidates), rate -4 ln 3 (3 candidates), coordination -4 ln 3 (3 unordered pairs)
capture scalar drop __nwdynam_sequential
mata:
xll = .
xg = .
XS = DynamState()
XS.init((1,2,1 \ 2,1,2 \ 1,2,3 \ 1,3,4), 3)
dynam_choice_loglik_grad_unit1((0,0,0), &XS, xll, xg)
assert(abs(xll - (-4*ln(2))) < 1e-12)
assert(max(abs(xg :- (0, 0.5, 0))) < 1e-12)
dynam_rate_loglik_grad_unit1((0,0), &XS, xll, xg)
assert(abs(xll - (-4*ln(3))) < 1e-12)
assert(max(abs(xg :- (4/3, 1/3))) < 1e-12)
// the Mata multi engines on the same sequence: the same numbers
xa = (1, 1, 1, J(1, 14, 0))
dynam_choice_loglik_grad_multi((0,0,0), xa, J(7,3,0), (1e300,1e300), (0,0,0,0), 0, J(0,0,0), J(0,0,0), &XS, xll, xg)
assert(abs(xll - (-4*ln(2))) < 1e-12)
assert(max(abs(xg :- (0, 0.5, 0))) < 1e-12)
dynam_rate_loglik_grad_multi((0,0), (1,1,0,0,0), J(1,3,0), J(1,3,0), 1e300, 1e300, 0, 0, 0, &XS, xll, xg)
assert(abs(xll - (-4*ln(3))) < 1e-12)
assert(max(abs(xg :- (4/3, 1/3))) < 1e-12)
// one effect alone: the 1 x 1 path
xa = J(1, 17, 0)
xa[2] = 1
dynam_choice_loglik_grad_multi((0), xa, J(7,3,0), (1e300,1e300), (0,0,0,0), 0, J(0,0,0), J(0,0,0), &XS, xll, xg)
assert(cols(xg) == 1)
assert(abs(xg - 0.5) < 1e-12)
// coordination: uniform over the 3 unordered pairs at theta = 0
dynam_coord_loglik_grad_multi((0), (1,0,0,0,0,0,0,0,0,0,0,0), J(7,3,0), J(0,0,0), &XS, xll, xg)
assert(abs(xll - (-4*ln(3))) < 1e-12)
// simultaneous events 2 and 3
XT = DynamState()
XT.init((1,2,1 \ 2,1,2 \ 1,2,2 \ 1,3,3), 3)
assert(XT.ngroups == 3)
dynam_choice_loglik_grad_unit1((0,0,0), &XT, xll, xg)
assert(max(abs(xg :- (0, 0, 0))) < 1e-12)
dynam_rate_loglik_grad_unit1((0,0), &XT, xll, xg)
assert(max(abs(xg :- (2/3, 2/3))) < 1e-12)
end
* option sequential (read by DynamState::init from the scalar the ado sets)
scalar __nwdynam_sequential = 1
mata:
XQ = DynamState()
XQ.init((1,2,1 \ 2,1,2 \ 1,2,2 \ 1,3,3), 3)
assert(XQ.ngroups == 4)
dynam_choice_loglik_grad_unit1((0,0,0), &XQ, xll, xg)
assert(max(abs(xg :- (0, 0.5, 0))) < 1e-12)
dynam_rate_loglik_grad_unit1((0,0), &XQ, xll, xg)
assert(max(abs(xg :- (4/3, 1/3))) < 1e-12)
end
scalar drop __nwdynam_sequential
di as text "  PASS: hand-derived scores, sequential and simultaneous"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: closed-form MLE on W40"
di as text "{hline 60}"
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
xt_cov w40
assert _N == 4
gen double x01 = (_n >= 3)
gen double x4 = (_n == 4)
nwdynam w40, submodel(rate) ego(x01) seed(1)
assert "`e(effects)'" == "ego"
xt_cf "rate ego(x01)" -55.0000715716 -0.302280871873 0.319846510504
nwdynam w40, submodel(rate) intercept ego(x01) seed(1)
assert "`: colnames e(b)'" == "Intercept ego"
xt_cf "rate intercept ego(x01)" -93.9873592523 -1.22121461076 0.208514414057 -0.302280871873 0.319846510504
nwdynam w40, alter(x4) seed(1)
assert "`e(effects)'" == "alter"
xt_cf "choice alter(x4)" -43.2266914797 -0.496436886314 0.431665847302
matrix XNX = (0,1,0,0 \ 0,0,1,0 \ 0,0,0,1 \ 1,0,0,0)
nwset, mat(XNX) directed name(nxt) nooutput
nwdynam w40, tie(nxt) seed(1)
assert "`e(effects)'" == "tie"
xt_cf "choice tie(nxt)" -40.4753676513 -1.04145387483 0.44280744277
* the same closed forms hold with tied times (the statistics are static): ceil(t/4) and option sequential
nwclear
xt_w40
quietly replace t = ceil(t/4)
nwset sender receiver, eventtime(t) directed name(wq)
xt_cov wq
gen double x01 = (_n >= 3)
gen double x4 = (_n == 4)
nwdynam wq, submodel(rate) ego(x01) seed(1)
xt_cf "rate ego(x01), ties" -55.0000715716 -0.302280871873 0.319846510504
nwdynam wq, submodel(rate) ego(x01) sequential seed(1)
xt_cf "rate ego(x01), ties, sequential" -55.0000715716 -0.302280871873 0.319846510504
nwdynam wq, alter(x4) seed(1)
xt_cf "choice alter(x4), ties" -43.2266914797 -0.496436886314 0.431665847302
* intercept with ties: the waiting time inside a time point is 0, so T is still t_last - t_first = 10 - 1 = 9:
* Intercept = ln(23/(2*9)), ego = ln(17/23), ll = 23 ln(23/18) + 17 ln(17/18) - 40
local i0 = ln(23/18)
local i1 = ln(17/23)
local llq = 23*ln(23/18) + 17*ln(17/18) - 40
nwdynam wq, submodel(rate) intercept ego(x01) seed(1)
xt_cf "rate intercept ego(x01), ties" `llq' `i0' `=sqrt(1/23)' `i1' `=sqrt(1/17 + 1/23)'

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: the smallest data"
di as text "{hline 60}"
* two events on 3 actors: rate ego(x) with x = (1,0,0), events 1->2, 2->3: a = b = 1, n1 = 1, n0 = 2:
* theta = ln 2, ll = ln(1/2) + ln(1/4), se = sqrt(2)
nwclear
clear
input long sender long receiver double t
1 2 1
2 3 2
end
nwset sender receiver, eventtime(t) directed name(two)
xt_cov two
assert _N == 3
gen double x = (_n == 1)
nwdynam two, submodel(rate) ego(x) seed(1)
assert e(N) == 2 & e(nodes) == 3
xt_cf "2 events, rate ego" -2.07944154168 0.69314718056 1.41421356237
* two events on 3 actors: choice alter(x) with x = (0,0,1), events 1->3, 2->1: one candidate with x = 1 among two in
* both events, a = b = 1: theta = 0, ll = 2 ln(1/2), se = sqrt(2)
nwclear
clear
input long sender long receiver double t
1 3 1
2 1 2
end
nwset sender receiver, eventtime(t) directed name(two2)
xt_cov two2
gen double x = (_n == 3)
nwdynam two2, alter(x) seed(1)
assert abs(_b[alter]) < 1e-4
assert reldif(e(ll), -1.38629436112) < 1e-8
assert reldif(_se[alter], 1.41421356237) < 1e-3
di as text "  PASS: 2 events, choice alter: b = 0, ll = 2 ln(1/2)"
* two actors: the choice of the receiver is forced at every event; the likelihood is 0 for every coefficient.
nwclear
clear
input long sender long receiver double t
1 2 1
2 1 2
1 2 3
end
nwset sender receiver, eventtime(t) directed name(pair)
capture noisily nwdynam pair, seed(1)
if _rc == 0 {
	assert abs(e(ll)) < 1e-12
	di as text "  PASS: 2 actors, choice: ll = 0 (unidentified: `e(unidentified)')"
}
else {
	assert _rc < 3000 | _rc > 3999
	di as text "  PASS: 2 actors, choice refused with r(" _rc ")"
}
* rate on 2 actors with ego(x), x = (1,0): a = 2 (actor 1 sends twice), b = 1, n1 = n0 = 1: theta = ln 2,
* ll = 2 ln(2/3) + ln(1/3), se = sqrt(3/2)
xt_cov pair
gen double x = (_n == 1)
nwdynam pair, submodel(rate) ego(x) seed(1)
xt_cf "2 actors, rate ego" -1.90954250488 0.69314718056 1.22474487139

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: many actors, few events (60 actors, 30 events; event i: i -> i + 30)"
di as text "{hline 60}"
* rate ego(x), x = 1 for actors 1..10: a = 10, b = 20, n1 = 10, n0 = 50: theta = ln 2.5,
* ll = 10 ln(10/300) + 20 ln(20/1500), se = sqrt(0.15). Actors 31..60 never send (they are in the risk set).
nwclear
clear
quietly set obs 30
gen long sender = _n
gen long receiver = _n + 30
gen double t = _n
nwset sender receiver, eventtime(t) directed name(big)
xt_cov big
assert _N == 60
gen double x = (_n <= 10)
nwdynam big, submodel(rate) ego(x) seed(1)
assert e(nodes) == 60
xt_cf "60 actors, rate ego" -120.361736087 0.916290731874 0.387298334621

di "ALL PASS"
