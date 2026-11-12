cscript

* xt_dyn_edges.do -- DRAFT, not yet run.
* Event-data edge grid for nwdynam: unsorted rows, time 0 and negative times, monotone time transforms (no window,
* no intercept: only the order enters), tied times (simultaneous as in goldfish; option sequential), shuffled rows
* inside a time point, all events at one time point, string / mixed / offset / permuted actor ids, valued event
* networks, two-mode label variants (numeric overlap, string overlap), undirected coordination events listed as
* (a,b) or (b,a). Reference values: goldfish 1.6.12 on the W40 toy sequence (R code in xt_NOTES.md), or invariances.
* Expected runtime: about 40 seconds (Stata 17 BE).
* Run from the repository root; the test runner sets the adopath (core, this repository, the plugin folder).
* No sysdir/adopath here.

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
	quietly gen long row = _n
end

capture program drop xt_gf
program xt_gf
	* xt_gf "label" ll b1 [b2 ...]: goldfish reference; coefficients to 2e-3, ll to 1e-5 relative (as test_nwdynam_ties.do)
	gettoken label 0 : 0
	gettoken ll 0 : 0
	local k = colsof(e(b))
	forvalues j = 1/`k' {
		gettoken ref 0 : 0
		assert abs(el(e(b), 1, `j') - (`ref')) < 2e-3 * (1 + abs(`ref'))
	}
	assert reldif(e(ll), `ll') < 1e-5
	di as text "  PASS: `label' = goldfish (ll " %12.6f e(ll) ")"
end

capture program drop xt_same
program xt_same
	* xt_same "label" tol prefix: e(b) and e(ll) equal the saved fit `prefix'B / `prefix'LL
	args label tol p
	assert colsof(e(b)) == colsof(`p'B)
	assert mreldif(e(b), `p'B) <= `tol'
	assert reldif(e(ll), scalar(`p'LL)) <= `tol'
	di as text "  PASS: `label'"
end

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: W40 reference fits (choice native engine, rate native engine, rate with intercept)"
di as text "{hline 60}"
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
nwdynam w40, seed(1)
assert "`e(effects)'" == "inertia recip indeg"
assert "`e(ties)'" == "none" & e(ntimes) == 40
xt_gf "choice default" -38.1300248 1.0497675 1.0747602 -0.2707421
matrix CB = e(b)
scalar CLL = e(ll)
nwdynam w40, submodel(rate) seed(1)
assert "`e(effects)'" == "indeg outdeg"
xt_gf "rate default" -55.41977877 0.02063502 0.09236512
matrix RB = e(b)
scalar RLL = e(ll)
nwdynam w40, submodel(rate) intercept seed(1)
assert "`: colnames e(b)'" == "Intercept indeg outdeg"
xt_gf "rate intercept" -94.37558055 -1.259479036 -0.067194365 0.009146213
matrix IB = e(b)
scalar ILL = e(ll)
* the Mata engine (any window forces it) agrees with the native engine: a window longer than the whole log
nwdynam w40, inertia recip indeg inertiawindow(1e9) seed(1)
xt_gf "choice, Mata engine" -38.1300248 1.0497675 1.0747602 -0.2707421
nwdynam w40, submodel(rate) indeg outdeg outdegwindow(1e9) seed(1)
xt_gf "rate, Mata engine" -55.41977877 0.02063502 0.09236512

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: unsorted rows, time 0, negative times, monotone time transforms"
di as text "{hline 60}"
local k 0
foreach texp in "t" "t - 20" "t - 1" "t * 1000" "t^3" "exp(t/10)" "-1/t" {
	local ++k
	nwclear
	xt_w40
	quietly replace t = `texp'
	set seed 4
	gen double u = runiform()
	sort u
	nwset sender receiver, eventtime(t) directed name(wt`k')
	nwdynam wt`k', seed(1)
	xt_same "choice, rows shuffled, time = `texp'" 1e-12 C
	nwdynam wt`k', submodel(rate) seed(1)
	xt_same "rate, rows shuffled, time = `texp'" 1e-12 R
}
* intercept: only the differences of the times enter; a shift leaves the fit unchanged
nwclear
xt_w40
quietly replace t = t + 1000
nwset sender receiver, eventtime(t) directed name(wshift)
nwdynam wshift, submodel(rate) intercept seed(1)
xt_same "rate intercept, times + 1000" 1e-10 I

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: tied times (ceil(t/4): 10 time points) - simultaneous, as goldfish"
di as text "{hline 60}"
nwclear
xt_w40
quietly replace t = ceil(t/4)
tempfile tq
quietly save `tq'
nwset sender receiver, eventtime(t) directed name(wq)
nwdynam wq, seed(1)
assert e(ntimes) == 10 & e(N) == 40
assert "`e(ties)'" == "simultaneous"
xt_gf "choice, ties" -38.9031366290 0.9853781955 1.0003945242 -0.3635616050
matrix TCB = e(b)
scalar TCLL = e(ll)
nwdynam wq, submodel(rate) seed(1)
xt_gf "rate, ties" -54.7027768468 -0.2559552691 -0.5005636724
matrix TRB = e(b)
scalar TRLL = e(ll)
nwdynam wq, submodel(rate) intercept seed(1)
xt_gf "rate intercept, ties" -34.19667944338 0.64072116773 -0.05028715403 -0.27224831861
* rows shuffled inside the time points: the same fit (all events of a time point see the same state)
nwclear
use `tq', clear
set seed 12
gen double u = runiform()
sort t u
nwset sender receiver, eventtime(t) directed name(wqs)
nwdynam wqs, seed(1)
xt_same "choice, ties, rows shuffled within time points" 1e-9 TC
nwdynam wqs, submodel(rate) seed(1)
xt_same "rate, ties, rows shuffled within time points" 1e-9 TR
* option sequential on tied data = the plain fit of distinct times in the listed order (here: W40 itself)
nwclear
use `tq', clear
sort row
nwset sender receiver, eventtime(t) directed name(wqq)
nwdynam wqq, sequential seed(1)
assert "`e(ties)'" == "sequential"
xt_same "sequential ties = distinct times" 1e-12 C
nwdynam wqq, submodel(rate) sequential seed(1)
xt_same "rate, sequential ties = distinct times" 1e-12 R
* sequential on data WITHOUT ties changes nothing
nwdynam w40, sequential seed(1)
assert "`e(ties)'" == "none"
xt_same "sequential without ties" 1e-12 C

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: all 40 events at one time point"
di as text "{hline 60}"
* every event sees the empty network: inertia, recip, indeg are 0 for every candidate, the likelihood is flat and
* equals the null value -40 ln 3 = -43.9444915467 (choice) and -40 ln 4 = -55.4517744448 (rate).
* Expected: a fit with that ll, or a clean refusal. (That every effect should then be in e(unidentified) is checked
* in xt_dyn_suspected.do, S8.)
nwclear
xt_w40
quietly replace t = 0
nwset sender receiver, eventtime(t) directed name(wzero)
foreach sm in choice rate {
	capture noisily nwdynam wzero, submodel(`sm') seed(1)
	local rc = _rc
	if `rc' == 0 {
		assert e(ntimes) == 1
		local null = cond("`sm'" == "choice", -43.9444915467, -55.4517744448)
		assert reldif(e(ll), `null') < 1e-10
		di as text "  PASS: `sm', one time point: ll = null ll (unidentified: `e(unidentified)')"
	}
	else {
		assert `rc' < 3000 | `rc' > 3999
		di as text "  PASS: `sm', one time point refused with r(`rc')"
	}
}

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block E: actor ids"
di as text "{hline 60}"
nwclear
xt_w40
gen str3 ss = "A" + string(sender)
gen str3 rs = "A" + string(receiver)
nwset ss rs, eventtime(t) directed name(wstr)
nwdynam wstr, seed(1)
xt_same "string ids" 1e-12 C
nwclear
xt_w40
gen str1 rs = string(receiver)
nwset sender rs, eventtime(t) directed name(wmix)
nwdynam wmix, seed(1)
xt_same "numeric sender, string receiver" 1e-12 C
local k 0
foreach map in "100 + @" "cond(@==1,1,cond(@==2,5,cond(@==3,17,1000)))" "10000000 + 2*@ - 1" "-10 + @" {
	local ++k
	nwclear
	xt_w40
	local es = subinstr("`map'", "@", "sender", .)
	local er = subinstr("`map'", "@", "receiver", .)
	quietly replace sender = `es'
	quietly replace receiver = `er'
	nwset sender receiver, eventtime(t) directed name(wid`k')
	nwdynam wid`k', seed(1)
	xt_same "ids `map'" 1e-12 C
}
* permuted labels (string ids that sort in reverse): estimates unchanged up to the order of floating-point sums
nwclear
xt_w40
gen str1 ss = word("d c b a", sender)
gen str1 rs = word("d c b a", receiver)
nwset ss rs, eventtime(t) directed name(wperm)
nwdynam wperm, seed(1)
xt_same "actor permutation, choice" 1e-6 C
nwdynam wperm, submodel(rate) seed(1)
xt_same "actor permutation, rate" 1e-6 R

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block F: a valued event network (third variable): unweighted effects ignore the values"
di as text "{hline 60}"
nwclear
xt_w40
gen double w = mod(_n, 3) + 1
nwset sender receiver w, eventtime(t) directed name(wval)
nwdynam wval, seed(1)
xt_same "valued event network, choice" 1e-12 C

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block G: two-mode label variants (dev/dynam_crosscheck_twomode_events.csv, 6 persons, 4 organisations)"
di as text "{hline 60}"
* goldfish reference (test_nwdynam_twomode_ado.do): choice inertia -0.73342, indeg 0.23018; rate outdeg 0.23797
nwclear
import delimited "dev/dynam_crosscheck_twomode_events.csv", clear varnames(1)
tempfile tm
quietly save `tm'
nwset sender receiver, twomode eventtime(time) name(tmstr)
nwdynam tmstr, seed(1)
assert e(nodes) == 10
assert abs(el(e(b),1,1) - (-0.73342)) < 0.001
assert abs(el(e(b),1,2) - 0.23018) < 0.001
matrix MB = e(b)
scalar MLL = e(ll)
local v 0
foreach variant in "numeric codes (overlapping ids)" "numeric codes, organisations + 100" "strings with overlapping labels" {
	local ++v
	nwclear
	use `tm', clear
	if `v' == 1 | `v' == 2 {
		encode sender, gen(sn)
		encode receiver, gen(rn)
		if `v' == 2 quietly replace rn = rn + 100
		nwset sn rn, twomode eventtime(time) name(tmv`v')
	}
	else {
		quietly replace sender = "X" + substr(sender, 2, .)
		quietly replace receiver = "X" + substr(receiver, 2, .)
		nwset sender receiver, twomode eventtime(time) name(tmv`v')
	}
	nwdynam tmv`v', seed(1)
	assert e(nodes) == 10
	xt_same "two-mode: `variant'" 1e-12 M
}

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block H: undirected coordination events, (a,b) or (b,a) listed"
di as text "{hline 60}"
* goldfish: choice_coordination inertia on W40 (undirected): 0.7290164305, ll -68.8462194594;
* with ceil(t/4): 0.4570586411, ll -70.5803895458
nwclear
xt_w40
nwset sender receiver, eventtime(t) undirected name(wu)
nwdynam wu, submodel(choice_coordination) inertia seed(1)
xt_gf "coordination inertia" -68.8462194594 0.7290164305
matrix UB = e(b)
scalar ULL = e(ll)
nwclear
xt_w40
gen long a = cond(mod(_n, 2), sender, receiver)
gen long b = cond(mod(_n, 2), receiver, sender)
nwset a b, eventtime(t) undirected name(wus)
nwdynam wus, submodel(choice_coordination) inertia seed(1)
xt_same "coordination, half of the pairs reversed" 1e-10 U
nwclear
xt_w40
quietly replace t = ceil(t/4)
nwset sender receiver, eventtime(t) undirected name(wuq)
nwdynam wuq, submodel(choice_coordination) inertia seed(1)
assert "`e(ties)'" == "simultaneous"
xt_gf "coordination inertia, ties" -70.5803895458 0.4570586411

di "ALL PASS"
