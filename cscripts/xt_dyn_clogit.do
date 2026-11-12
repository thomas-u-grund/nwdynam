cscript

* xt_dyn_clogit.do -- DRAFT, not yet run.
* Cross-route check: the choice sub-model is a conditional logit with one stratum per event (the sender's n - 1
* candidate receivers), the rate sub-model without intercept one with all n actors as candidates. This file computes
* the statistics with its own Mata code (written from the help file's definitions, not from unw_dynam.do), fits
* -clogit- and compares with nwdynam: structural, weighted, windowed and attribute effects, tied times
* (simultaneous: all events of a time point see the state before it) and option sequential, opportunities()
* (candidate rows dropped), and a permutation of the actor labels.
* Data: 6 actors, 200 events, generated below with rseed(20261113).
* Expected runtime: about 90 seconds (Stata 17 BE).
* Run from the repository root; the test runner sets the adopath (core, this repository, the plugin folder).
* No sysdir/adopath here.

capture program drop nwinstall
local nwcore : env NWCORE
if `"`nwcore'"' != "" quietly do `"`nwcore'/unw_core.do"'
quietly do unw_dynam.do

capture mata: mata drop xt_choice_long()
capture mata: mata drop xt_rate_long()
mata:
mata set matastrict off
// Choice: one row per event i and candidate j != sender. Events in the row order of E; tt = times (equal
// consecutive times form a group: every event of the group sees the state before the group). W = window.
// cx, cz = per-actor covariates (column vectors).
// columns: 1 ev 2 cand 3 chosen 4 inertia 5 recip 6 indeg 7 outdeg 8 trans 9 cycle 10 csender 11 creceiver
//          12 winertia 13 wrecip 14 windeg 15 woutdeg 16 nodetrans 17 inertiaw 18 recipw
//          19 same(cx) 20 diff(cz) 21 alter(cz) 22 egoalterint(cx,cz) 23 tertius(cz) 24 sim(cz)
real matrix xt_choice_long(real matrix E, real colvector tt, real scalar n, real scalar W, real colvector cx, real colvector cz)
{
	N = rows(E)
	A = J(n, n, 0)
	C = J(n, n, 0)
	LT = J(n, n, .)
	L = J(N * (n - 1), 24, .)
	row = 0
	g0 = 1
	for (i=1; i<=N; i++) {
		if (i > 1) {
			if (tt[i] != tt[i-1]) g0 = i
		}
		s = E[i,1]
		r = E[i,2]
		AA = A * A
		AtA = A' * A
		AAt = A * A'
		cs = colsum(A)
		rs = rowsum(A)'
		wcs = colsum(C)
		wrs = rowsum(C)'
		nt = rowsum(A :* AA)'
		for (j=1; j<=n; j++) {
			if (j == s) continue
			row = row + 1
			iw = 0
			if (LT[s,j] < .) iw = ((tt[i] - LT[s,j]) <= W)
			rw = 0
			if (LT[j,s] < .) rw = ((tt[i] - LT[j,s]) <= W)
			ter = 0
			if (cs[j] > 0) ter = (cz' * A[.,j]) / cs[j]
			L[row,.] = (i, j, (j == r), A[s,j], A[j,s], cs[j], rs[j], AA[s,j], AA[j,s], AtA[s,j], AAt[s,j],
				C[s,j], C[j,s], wcs[j], wrs[j], nt[j], iw, rw,
				(cx[j] == cx[s]), abs(cz[j] - cz[s]), cz[j], cx[s] * cz[j], ter, -abs(cz[j] - cz[s]))
		}
		lastg = (i == N)
		if (!lastg) lastg = (tt[i+1] != tt[i])
		if (lastg) {
			for (k=g0; k<=i; k++) {
				A[E[k,1], E[k,2]] = 1
				C[E[k,1], E[k,2]] = C[E[k,1], E[k,2]] + 1
				LT[E[k,1], E[k,2]] = tt[i]
			}
		}
	}
	return(L)
}
// Rate (no intercept): one row per event i and actor j (all n). columns: 1 ev 2 cand 3 chosen 4 indeg 5 outdeg
//          6 nodetrans 7 windeg 8 woutdeg 9 outdegw 10 indegw 11 ego(cz) 12 tertius(cz)
real matrix xt_rate_long(real matrix E, real colvector tt, real scalar n, real scalar W, real colvector cz)
{
	N = rows(E)
	A = J(n, n, 0)
	C = J(n, n, 0)
	LT = J(n, n, .)
	L = J(N * n, 12, .)
	row = 0
	g0 = 1
	for (i=1; i<=N; i++) {
		if (i > 1) {
			if (tt[i] != tt[i-1]) g0 = i
		}
		s = E[i,1]
		cs = colsum(A)
		rs = rowsum(A)'
		nt = rowsum(A :* (A * A))'
		WIN = (LT :< .) :& ((tt[i] :- LT) :<= W)
		wcs = colsum(C)
		wrs = rowsum(C)
		win_out = rowsum(WIN)
		win_in = colsum(WIN)
		for (j=1; j<=n; j++) {
			row = row + 1
			ter = 0
			if (cs[j] > 0) ter = (cz' * A[.,j]) / cs[j]
			L[row,.] = (i, j, (j == s), cs[j], rs[j], nt[j], wcs[j], wrs[j], win_out[j], win_in[j], cz[j], ter)
		}
		lastg = (i == N)
		if (!lastg) lastg = (tt[i+1] != tt[i])
		if (lastg) {
			for (k=g0; k<=i; k++) {
				A[E[k,1], E[k,2]] = 1
				C[E[k,1], E[k,2]] = C[E[k,1], E[k,2]] + 1
				LT[E[k,1], E[k,2]] = tt[i]
			}
		}
	}
	return(L)
}
end

* ---------------------------------------------------------------------------------------------------------------------
* data: 6 actors, 200 events with repetition and reciprocation tendencies; covariates cx (groups), cz (continuous)
mata:
rseed(20261113)
xn = 6
xN = 200
XE = J(xN, 2, 0)
lastr = J(1, xn, 0)
lasts = J(1, xn, 0)
for (i=1; i<=xN; i++) {
	s = runiformint(1, 1, 1, xn)
	u = runiform(1, 1)
	if (lastr[s] > 0 & u < 0.3) r = lastr[s]
	else if (lasts[s] > 0 & u < 0.5) r = lasts[s]
	else {
		r = runiformint(1, 1, 1, xn - 1)
		if (r >= s) r = r + 1
	}
	XE[i,.] = (s, r)
	lastr[s] = r
	lasts[r] = s
}
xcx = (1 \ 1 \ 2 \ 2 \ 3 \ 3)
xcz = (0.2 \ -1 \ 1.4 \ 0.5 \ -0.3 \ 2)
xt1 = (1::xN)
xt4 = ceil((1::xN) :/ 4)
end

capture program drop xt_events
program xt_events
	* the generated events; `1' = time expression (default _n)
	args texp
	if "`texp'" == "" local texp "_n"
	clear
	mata: st_local("nobs", strofreal(rows(XE)))
	quietly set obs `nobs'
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store(., ("sender", "receiver"), XE)
	quietly gen double t = `texp'
end

capture program drop xt_covs
program xt_covs
	args net
	nwload `net', xvars
	assert _N == 6
	quietly gen double cx = .
	quietly gen double cz = .
	mata: st_store(., ("cx", "cz"), (xcx, xcz))
end

capture program drop xt_keep
program xt_keep
	* save the current nwdynam fit as matrix `1'B and scalar `1'LL
	args p
	matrix `p'B = e(b)
	scalar `p'LL = e(ll)
end

capture program drop xt_cmp
program xt_cmp
	* xt_cmp prefix "nwdynam names" "clogit variables": saved fit against the clogit in memory
	args p nwn clv
	local k : word count `nwn'
	assert colsof(`p'B) == `k'
	forvalues q = 1/`k' {
		local a : word `q' of `nwn'
		local c : word `q' of `clv'
		local bn = el(`p'B, 1, colnumb(`p'B, "`a'"))
		assert abs(`bn' - _b[`c']) <= 1e-3 * (1 + abs(_b[`c']))
	}
	assert abs(scalar(`p'LL) - e(ll)) < 1e-5
	di as text "  PASS: `p': nwdynam (`nwn') = clogit (`clv'), ll " %12.6f e(ll)
end

capture program drop xt_long
program xt_long
	* long data in memory: `1' = choice or rate, `2' = Mata time vector, `3' = window
	args kind tv W
	clear
	if "`kind'" == "choice" {
		mata: XL = xt_choice_long(XE, `tv', xn, `W', xcx, xcz)
		mata: (void) st_addvar("double", ("ev","cand","chosen","inertia","recip","indeg","outdeg","trans","cycle","csender","creceiver","winertia","wrecip","windeg","woutdeg","nodetrans","inertiaw","recipw","same","diff","alter","egoalterint","tertius","sim"))
	}
	else {
		mata: XL = xt_rate_long(XE, `tv', xn, `W', xcz)
		mata: (void) st_addvar("double", ("ev","cand","chosen","indeg","outdeg","nodetrans","windeg","woutdeg","outdegw","indegw","ego","tertius"))
	}
	mata: st_addobs(rows(XL))
	mata: st_store(., ., XL)
	quietly count if chosen == 1
	assert r(N) == 200
end

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: nwdynam fits, distinct times"
di as text "{hline 60}"
nwclear
xt_events
nwset sender receiver, eventtime(t) directed name(gen)
xt_covs gen
nwdynam gen, seed(1)
xt_keep C1
nwdynam gen, outdeg trans cycle seed(1)
xt_keep C2
nwdynam gen, inertia commonsender commonreceiver nodetrans seed(1)
xt_keep C3
nwdynam gen, weightedinertia weightedrecip weightedindeg weightedoutdeg seed(1)
xt_keep C4
nwdynam gen, inertiawindow(3) recipwindow(3) indeg seed(1)
xt_keep C5
nwdynam gen, same(cx) diff(cz) alter(cz) seed(1)
xt_keep C6
nwdynam gen, inertia egoalterint(cx cz) tertius(cz) seed(1)
xt_keep C7
nwdynam gen, recip sim(cz) seed(1)
xt_keep C8
nwdynam gen, submodel(rate) seed(1)
xt_keep R1
nwdynam gen, submodel(rate) nodetrans weightedindeg weightedoutdeg seed(1)
xt_keep R2
nwdynam gen, submodel(rate) outdegwindow(3) indegwindow(3) seed(1)
xt_keep R3
nwdynam gen, submodel(rate) outdeg ego(cz) tertius(cz) seed(1)
xt_keep R4

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: opportunities() on the distinct-time data"
di as text "{hline 60}"
* candidate j of event i is unavailable when mod(i + 2*j, 5) == 0, unless j is the realised receiver
clear
mata: st_local("nobs", strofreal(rows(XE)))
quietly set obs `=`nobs' * 6'
gen long oev = ceil(_n / 6)
gen long oact = _n - 6 * (oev - 1)
gen long rec = .
mata: st_store(., "rec", XE[st_data(., "oev"), 2])
drop if mod(oev + 2*oact, 5) == 0 & oact != rec
nwdynam gen, inertia indeg opportunities(oev oact) seed(1)
xt_keep O1

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: nwdynam fits, tied times (ceil(i/4)), simultaneous and sequential"
di as text "{hline 60}"
nwclear
xt_events "ceil(_n/4)"
nwset sender receiver, eventtime(t) directed name(gq)
xt_covs gq
nwdynam gq, seed(1)
assert "`e(ties)'" == "simultaneous" & e(ntimes) == 50
xt_keep T1
nwdynam gq, outdeg trans cycle seed(1)
xt_keep T2
nwdynam gq, weightedinertia weightedrecip seed(1)
xt_keep T3
nwdynam gq, inertiawindow(2) recip seed(1)
xt_keep T4
nwdynam gq, submodel(rate) seed(1)
xt_keep T5
nwdynam gq, submodel(rate) indegwindow(2) nodetrans seed(1)
xt_keep T6
nwdynam gq, sequential seed(1)
assert "`e(ties)'" == "sequential"
assert mreldif(e(b), C1B) < 1e-12
nwdynam gq, submodel(rate) sequential seed(1)
assert mreldif(e(b), R1B) < 1e-12
di as text "  PASS: sequential ties = distinct times"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: clogit, distinct times"
di as text "{hline 60}"
xt_long choice xt1 3
quietly clogit chosen inertia recip indeg, group(ev)
xt_cmp C1 "inertia recip indeg" "inertia recip indeg"
quietly clogit chosen outdeg trans cycle, group(ev)
xt_cmp C2 "outdeg trans cycle" "outdeg trans cycle"
quietly clogit chosen inertia csender creceiver nodetrans, group(ev)
xt_cmp C3 "inertia commonsender commonreceiver nodetrans" "inertia csender creceiver nodetrans"
quietly clogit chosen winertia wrecip windeg woutdeg, group(ev)
xt_cmp C4 "inertia recip indeg outdeg" "winertia wrecip windeg woutdeg"
quietly clogit chosen inertiaw recipw indeg, group(ev)
xt_cmp C5 "inertia recip indeg" "inertiaw recipw indeg"
quietly clogit chosen same diff alter, group(ev)
xt_cmp C6 "same diff alter" "same diff alter"
quietly clogit chosen inertia egoalterint tertius, group(ev)
xt_cmp C7 "inertia egoalterint tertius" "inertia egoalterint tertius"
quietly clogit chosen recip sim, group(ev)
xt_cmp C8 "recip sim" "recip sim"
* opportunities: drop the unavailable candidate rows
quietly bysort ev: egen rec = max(cond(chosen == 1, cand, .))
drop if mod(ev + 2*cand, 5) == 0 & cand != rec
quietly clogit chosen inertia indeg, group(ev)
xt_cmp O1 "inertia indeg" "inertia indeg"
xt_long rate xt1 3
quietly clogit chosen indeg outdeg, group(ev)
xt_cmp R1 "indeg outdeg" "indeg outdeg"
quietly clogit chosen windeg woutdeg nodetrans, group(ev)
xt_cmp R2 "indeg outdeg nodetrans" "windeg woutdeg nodetrans"
quietly clogit chosen indegw outdegw, group(ev)
xt_cmp R3 "indeg outdeg" "indegw outdegw"
quietly clogit chosen outdeg ego tertius, group(ev)
xt_cmp R4 "outdeg ego tertius" "outdeg ego tertius"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block E: clogit, tied times (simultaneous)"
di as text "{hline 60}"
xt_long choice xt4 2
quietly clogit chosen inertia recip indeg, group(ev)
xt_cmp T1 "inertia recip indeg" "inertia recip indeg"
quietly clogit chosen outdeg trans cycle, group(ev)
xt_cmp T2 "outdeg trans cycle" "outdeg trans cycle"
quietly clogit chosen winertia wrecip, group(ev)
xt_cmp T3 "inertia recip" "winertia wrecip"
quietly clogit chosen inertiaw recip, group(ev)
xt_cmp T4 "inertia recip" "inertiaw recip"
xt_long rate xt4 2
quietly clogit chosen indeg outdeg, group(ev)
xt_cmp T5 "indeg outdeg" "indeg outdeg"
quietly clogit chosen indegw nodetrans, group(ev)
xt_cmp T6 "indeg nodetrans" "indegw nodetrans"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block F: permuted actor labels (1..6 -> 4,6,1,5,2,3) with permuted covariates"
di as text "{hline 60}"
mata:
xpi = (4 \ 6 \ 1 \ 5 \ 2 \ 3)
xpinv = invorder(xpi)
XEp = (xpi[XE[.,1]], xpi[XE[.,2]])
end
nwclear
clear
quietly set obs 200
gen long sender = .
gen long receiver = .
mata: st_store(., ("sender", "receiver"), XEp)
gen double t = _n
nwset sender receiver, eventtime(t) directed name(gp)
nwload gp, xvars
gen double cx = .
gen double cz = .
mata: st_store(., ("cx", "cz"), (xcx[xpinv], xcz[xpinv]))
nwdynam gp, seed(1)
assert mreldif(e(b), C1B) < 1e-6
assert reldif(e(ll), scalar(C1LL)) < 1e-9
nwdynam gp, same(cx) diff(cz) alter(cz) seed(1)
assert mreldif(e(b), C6B) < 1e-5
assert reldif(e(ll), scalar(C6LL)) < 1e-9
nwdynam gp, submodel(rate) outdeg ego(cz) tertius(cz) seed(1)
assert mreldif(e(b), R4B) < 1e-5
assert reldif(e(ll), scalar(R4LL)) < 1e-9
di as text "  PASS: label permutation leaves estimates unchanged"

di "ALL PASS"
