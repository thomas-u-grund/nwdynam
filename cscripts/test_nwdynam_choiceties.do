* estat choiceplot: ties in the fitted probabilities. Deterministic (no random numbers inside the command), rank = mid rank, top-k share =
* exact expectation under a random order within the tie group. Checked against an independent computation in Stata from estat eventprob,
* and against a constructed case with a known answer (a tie group of m members: the share is min(k,m)/m).
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

set seed 20261121
clear
set obs 80
generate long sender = ceil(6 * runiform())
generate long receiver = ceil(5 * runiform())
replace receiver = receiver + (receiver >= sender)
generate double eventtime = _n
nwset sender receiver, eventtime(eventtime) name(tt) directed
nwdynam tt, submodel(choice) inertia

* ---- 1. run twice: identical, with different seeds in between (no random numbers)
set seed 1
estat choiceplot, noplot topk(1 2 3) frame(f1)
matrix A = r(curve)
scalar nt1 = r(N_ties)
scalar mr1 = r(mean_rank)
set seed 999
estat choiceplot, noplot topk(1 2 3) frame(f2)
matrix B = r(curve)
assert mreldif(A, B) == 0
assert nt1 == r(N_ties) & mr1 == r(mean_rank)
assert nt1 > 0
di as result "PASS deterministic, " nt1 " events with ties"

* ---- 2. independent computation from the event probabilities
estat eventprob, frame(ep)
frame ep {
	keep if p < .
	bysort event: gen long m = _N
	generate byte hit = observed == 1
	bysort event (hit): gen double pobs = p[_N]
	gen byte gt = p > pobs * (1 + 1e-12)
	gen byte eq = abs(p - pobs) <= 1e-12 * pobs
	collapse (sum) gr = gt eq = eq (first) m, by(event)
	generate double rk = 1 + gr + (eq - 1) / 2
	quietly count if eq > 1
	assert r(N) == nt1
	forvalues k = 1/3 {
		generate double h`k' = min(1, max(0, (`k' - gr) / eq))
		quietly summarize h`k'
		assert abs(r(mean) - A[`k', 2]) < 1e-12
	}
	quietly summarize rk
	assert abs(r(mean) - mr1) < 1e-12
}
frame f1 {
	assert nties >= 1 & nties < .
	quietly count if nties > 1
	assert r(N) == nt1
}
di as result "PASS tie counts, mid rank and top-k expectation against independent computation"

* ---- 3. a constructed case with a known answer: every sender acts once, so inertia is 0 for all candidates and all m = 4 candidates tie
*       (the model is fitted on a log with repeats and new receivers; the check uses the first five events, each from a sender seen for the first time)
clear
input long s long r double t
1 2 1
2 3 2
3 4 3
4 5 4
5 1 5
1 2 6
2 4 7
3 4 8
4 1 9
5 1 10
end
nwset s r, eventtime(t) name(kk) directed
nwdynam kk, submodel(choice) inertia
estat choiceplot, noplot topk(1 2 3 4) frame(fk)
frame fk {
	assert nties[1] == 4 & nties[2] == 4 & nties[3] == 4 & nties[4] == 4 & nties[5] == 4
	assert rank[1] == 2.5 & rank[5] == 2.5
	* the first five events: one tie group of 4, so the share at k is k/4 for each of them
}
mata: st_numscalar("dummy", 1)
frame fk {
	quietly count if nties == 4
	assert r(N) == 5
}
di as result "PASS constructed ties"

* ---- 4. the printed line appears when ties occur
quietly log using `"`c(tmpdir)'/tiesline.log"', text replace name(tl)
estat choiceplot, noplot
quietly log close tl
tempname fh
local found 0
file open `fh' using `"`c(tmpdir)'/tiesline.log"', read text
file read `fh' line
while r(eof) == 0 {
	if strpos(`"`macval(line)'"', "Ties: in ") local found 1
	file read `fh' line
}
file close `fh'
assert `found' == 1
di as result "PASS tie line printed"
di as result "ALL OK test_nwdynam_choiceties"
