cscript

* estat choiceplot and estat rateplot (extensions: goldfish has no such plots): numbers against an independent enumeration in Mata
* (choice with inertia, rate with indeg, rate with intercept and indeg, on a small random event log with distinct times; the
* probabilities are written out from the coefficients, not taken from the plugin), the graphs exist and can be exported, the
* choice_coordination case (risk set = all pairs, identities), and the error paths.
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

mata:
mata set matastrict off
// events s r t (sorted by t, distinct times), n actors, b the coefficient vector; type 1 = choice/inertia, 2 = rate/indeg (offset 0 = no intercept coefficient)
void enum_test(real matrix ev, real scalar n, real scalar type, real rowvector b, real scalar K) {
	real scalar i, j, s, r, m, gr, eq, pobs, k
	real matrix cnt
	real colvector x, p, ok, srt, rk, mm
	real matrix curve, act
	cnt = J(n, n, 0)
	rk = J(rows(ev), 1, .); mm = rk
	curve = J(K, 3, 0)
	act = J(n, 1, 0)
	for (i = 1; i <= rows(ev); i++) {
		s = ev[i, 1]; r = ev[i, 2]
		if (type == 1) {
			x = (cnt[s, .]' :> 0)
			ok = J(n, 1, 1); ok[s] = 0
			p = exp(b[1] :* x) :* ok
			p = p / sum(p)
			pobs = p[r]
			p = select(p, ok)
			m = rows(p)
			gr = sum(p :> pobs * (1 + 1e-12)); eq = sum(abs(p :- pobs) :<= 1e-12 * pobs)
			rk[i] = 1 + gr + (eq - 1) / 2; mm[i] = m
			srt = runningsum(sort(p, -1))
			for (k = 1; k <= K; k++) {
				curve[k, 1] = curve[k, 1] + min((1, max((0, (k - gr) / eq))))
				curve[k, 2] = curve[k, 2] + min((k, m)) / m
				curve[k, 3] = curve[k, 3] + srt[min((k, m))]
			}
		}
		else {
			x = colsum(cnt :> 0)'
			p = exp(b[cols(b)] :* x)
			p = p / sum(p)
			act = act + p
		}
		cnt[s, r] = cnt[s, r] + 1
	}
	if (type == 1) {
		st_matrix("E_curve", curve / rows(ev))
		st_numscalar("E_mean", mean(rk))
		st_numscalar("E_median", mm[1] < . ? (sort(rk, 1)[ceil(rows(rk) / 2)] + sort(rk, 1)[floor(rows(rk) / 2) + 1]) / 2 : .)
		st_numscalar("E_uniform", mean((mm :+ 1) / 2))
	}
	else st_matrix("E_act", act)
}
end

set seed 20261120
clear
set obs 60
generate long sender = ceil(6 * runiform())
generate long receiver = ceil(5 * runiform())
replace receiver = receiver + (receiver >= sender)
generate double eventtime = _n + runiform() / 4
tempfile evfile
sort eventtime
save `evfile'
nwset sender receiver, eventtime(eventtime) name(toy) directed
use `evfile', clear
sort eventtime
mata: EV = st_data(., "sender receiver")

* ---- 1. choice sub-model: rank of the realised receiver
nwdynam toy, submodel(choice) inertia
matrix b = e(b)
mata: enum_test(EV, 6, 1, st_matrix("b"), 20)
estat choiceplot, topk(1 3 5) name(gchoice) frame(rkf) title("Passed title") xlabel(1(1)5) scheme(s2color)
assert abs(r(mean_rank) - E_mean) < 1e-10
assert abs(r(median_rank) - E_median) < 1e-10
assert abs(r(mean_rank_uniform) - E_uniform) < 1e-10
assert r(mean_riskset) == 5
assert r(N) == e(N)
matrix C = r(curve)
assert rowsof(C) == 20
forvalues k = 1/20 {
	assert abs(C[`k', 2] - E_curve[`k', 1]) < 1e-10
	assert abs(C[`k', 3] - E_curve[`k', 2]) < 1e-10
	assert abs(C[`k', 4] - E_curve[`k', 3]) < 1e-10
}
matrix T = r(topk)
assert rowsof(T) == 3
assert T[1, 1] == 1 & T[2, 1] == 3 & T[3, 1] == 5
assert abs(T[2, 2] - E_curve[3, 1]) < 1e-10
assert abs(T[3, 3] - 1) < 1e-12
* the curves are non-decreasing and reach 1 at the size of the risk set
forvalues k = 2/20 {
	assert C[`k', 2] >= C[`k' - 1, 2] - 1e-12
}
assert abs(C[5, 2] - 1) < 1e-12 & abs(C[5, 4] - 1) < 1e-12 & abs(C[5, 3] - 1) < 1e-12
* the per-event frame
frame rkf {
	assert _N == 60
	assert m == 5
	assert inrange(rank, 1, 5)
	assert pobs > 0 & pobs < 1
	quietly summarize rank
	assert abs(r(mean) - E_mean) < 1e-10
}
* graph exists and exports
graph describe gchoice
tempfile png
graph export `"`png'"', as(png) name(gchoice) width(600) replace
confirm file `"`png'"'
* default top-k table
estat choiceplot, noplot
matrix T = r(topk)
assert rowsof(T) == 4 & T[4, 1] == 10
* maxrank() larger than the risk set: curves stay at 1
estat choiceplot, maxrank(8) noplot
assert rowsof(r(curve)) == 10
di as result "PASS choice plot numbers agree with the enumeration"

* ---- 2. rate sub-model without intercept: expected number of events per actor
nwdynam toy, submodel(rate) indeg
matrix b = e(b)
mata: enum_test(EV, 6, 2, st_matrix("b"), 1)
estat rateplot, name(grate) frame(raf) label(2) title("Passed title") xtitle("Passed x")
assert r(N) == 60
assert r(N_actors) == 6
assert "`r(over)'" == "actor"
frame raf {
	assert _N == 6
	quietly summarize observed
	assert r(sum) == 60
	quietly summarize expected
	assert abs(r(sum) - 60) < 1e-9
	sort actorid
	forvalues a = 1/6 {
		assert abs(expected[`a'] - E_act[`a', 1]) < 1e-9
	}
	* observed = number of events as sender
}
forvalues a = 1/6 {
	quietly count if sender == `a'
	frame raf: assert observed[`a'] == r(N)
}
graph describe grate
estat rateplot, over(time) periods(4) topactors(3) name(gtime) frame(rpf) ylabel(, format(%9.1f)) xtitle("Passed x")
graph describe gtime
frame rpf {
	assert _N == 6 * 4
	quietly summarize observed
	assert r(sum) == 60
	quietly summarize expected
	assert abs(r(sum) - 60) < 1e-9
	assert tlo < thi
	collapse (sum) observed expected, by(period)
	assert _N == 4
}
* the periods split the events as counted directly
quietly summarize eventtime
local t0 = r(min)
local t1 = r(max)
forvalues q = 1/4 {
	quietly count if min(4, 1 + floor((eventtime - `t0') / (`t1' - `t0') * 4)) == `q'
	frame rpf: assert observed[`q'] == r(N)
}
di as result "PASS rate plot numbers agree with the enumeration (no intercept)"

* ---- 3. rate sub-model with intercept: share of the total hazard
nwdynam toy, submodel(rate) intercept indeg
matrix b = e(b)
mata: enum_test(EV, 6, 2, st_matrix("b"), 1)
estat rateplot, noplot frame(ra2)
frame ra2 {
	sort actorid
	forvalues a = 1/6 {
		assert abs(expected[`a'] - E_act[`a', 1]) < 1e-9
	}
}
estat rateplot, over(time) periods(3) name(gtime2) frame(rp2)
graph describe gtime2
di as result "PASS rate plot with intercept"

* ---- 4. choice_coordination: the risk set is the set of pairs
clear
set obs 50
generate long sender = ceil(5 * runiform())
generate long receiver = ceil(4 * runiform())
replace receiver = receiver + (receiver >= sender)
generate double eventtime = _n + runiform() / 4
nwset sender receiver, eventtime(eventtime) name(toyu) undirected
nwdynam toyu, submodel(choice_coordination) inertia
estat choiceplot, topk(1 2) frame(rku) name(gcoord)
assert r(mean_riskset) == 10
assert r(mean_rank) >= 1 & r(mean_rank) <= 10
assert abs(r(mean_rank_uniform) - 5.5) < 1e-12
graph describe gcoord
frame rku {
	assert _N == 50
	assert m == 10
	assert inrange(rank, 1, 10)
}
di as result "PASS choice_coordination rank plot"

* ---- 5. error paths
nwdynam toy, submodel(rate) indeg
capture estat choiceplot
assert _rc == 198
nwdynam toy, submodel(choice) inertia
capture estat rateplot
assert _rc == 198
capture estat choiceplot, maxrank(0)
assert _rc == 198
capture estat choiceplot, topk(0)
assert _rc != 0
capture estat choiceplot, frame(rkf)
assert _rc == 110
estat choiceplot, frame(rkf) replace noplot
capture estat choiceplot, maxrows(10)
assert _rc != 0
nwdynam toy, submodel(rate) indeg
capture estat rateplot, over(foo)
assert _rc == 198
capture estat rateplot, periods(0)
assert _rc == 198
capture estat rateplot, topactors(0)
assert _rc == 198
capture estat rateplot, frame(raf)
assert _rc == 110
* the fit in memory must be the one of the estimates
estimates store old
nwdynam toy, submodel(rate) indeg outdeg
estimates restore old
capture estat rateplot, noplot
assert _rc == 498
capture estat choiceplot, noplot
assert _rc == 498 | _rc == 198
ereturn clear
capture estat rateplot
assert _rc == 301
capture estat choiceplot
assert _rc == 301
di as result "PASS error paths"
di as result "ALL TESTS PASSED (choiceplot, rateplot)"
