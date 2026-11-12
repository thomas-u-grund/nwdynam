cscript

* xt_dyn_estimation.do -- DRAFT, not yet run.
* Determinism, stored results and postestimation for nwdynam:
*   seed reproducibility (Mata multi engine: random starts) and seed independence of the optimum; native unit1
*   engines (default effect sets) do not depend on the seed; every documented e() entry (help file, Stored results);
*   nwdynam version (r(reference) goldfish 1.6.12, r(license) GPL-3, r(version) = e(version)); ereturn display,
*   estat vce, lincom/test/nlcom identities, estimates store/restore/replay on every result shape (choice native,
*   choice multi, rate, rate with intercept, choice_coordination, two-mode); the data in memory, c(frame), the
*   event list and the scalar/matrix name space are unchanged by a fit.
* References: goldfish 1.6.12 on the W40 sequence (R code in xt_NOTES.md).
* Expected runtime: about 60 seconds (Stata 17 BE).
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
end

capture program drop xt_post
program xt_post
	* postestimation identities on the result in memory; `1' = label
	args label
	tempname b V Vt
	matrix `b' = e(b)
	matrix `V' = e(V)
	local k = colsof(`b')
	local names : colnames `b'
	assert "`names'" == "`e(effects)'" | "`names'" == "Intercept `e(effects)'"
	assert rowsof(`V') == `k' & colsof(`V') == `k'
	matrix `Vt' = `V''
	assert mreldif(`V', `Vt') < 1e-10
	forvalues j = 1/`k' {
		assert el(`V', `j', `j') > 0
	}
	ereturn display
	estat vce
	local n1 : word 1 of `names'
	lincom `n1'
	assert abs(r(estimate) - el(`b', 1, 1)) < 1e-12
	assert abs(r(se) - sqrt(el(`V', 1, 1))) < 1e-12
	test `n1'
	assert abs(r(chi2) - el(`b', 1, 1)^2 / el(`V', 1, 1)) < 1e-8 * (1 + r(chi2))
	assert r(df) == 1
	if `k' >= 2 {
		local n2 : word 2 of `names'
		lincom `n1' - `n2'
		assert abs(r(estimate) - (el(`b',1,1) - el(`b',1,2))) < 1e-12
		assert abs(r(se) - sqrt(el(`V',1,1) + el(`V',2,2) - 2*el(`V',1,2))) < 1e-10
		test `n1' `n2'
		assert r(df) == 2
		nlcom exp(_b[`n1'])
		assert abs(el(r(b),1,1) - exp(el(`b',1,1))) < 1e-10 * (1 + exp(el(`b',1,1)))
	}
	estimates store xt_st
	ereturn clear
	estimates restore xt_st
	assert mreldif(e(b), `b') == 0
	assert "`e(cmd)'" == "nwdynam"
	// POLICY N1 (owner: open): nwdynam has no replay; `nwdynam` alone refits the current network, so estimates replay
	// fails (or refits). Recorded, not counted.
	capture noisily estimates replay xt_st
	if _rc di as text "  NOTE N1: estimates replay gives r(`=_rc') (no replay in nwdynam)"
	estimates drop xt_st
	di as text "  PASS: postestimation, `label'"
end

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: seed - native engines do not use it, the multi engine reproduces with it"
di as text "{hline 60}"
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
* default choice and rate sets run on the native unit1 engines with fixed starts: identical for any seed
nwdynam w40, seed(1)
matrix C1 = e(b)
scalar CLL = e(ll)
foreach s in 2 99 123456 {
	nwdynam w40, seed(`s')
	assert mreldif(e(b), C1) == 0 & e(ll) == scalar(CLL)
}
set seed 42
nwdynam w40
assert mreldif(e(b), C1) == 0
nwdynam w40, submodel(rate) seed(1)
matrix R1 = e(b)
nwdynam w40, submodel(rate) seed(7)
assert mreldif(e(b), R1) == 0
di as text "  PASS: native engines independent of seed()"
* 8-effect choice model (multi engine): goldfish 2.3222278923 2.6066294426 -1.8832365523 -1.3145112385
* 1.2728299833 1.1697045224 0.3808007836 -0.6188382957, ll -33.1426170323 (effects listed in goldfish order; nwdynam
* stores inertia recip indeg outdeg trans cycle commonsender commonreceiver in the same order)
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(11)
assert "`e(effects)'" == "inertia recip indeg outdeg trans cycle commonsender commonreceiver"
matrix M11 = e(b)
scalar MLL = e(ll)
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(11)
assert mreldif(e(b), M11) == 0 & e(ll) == scalar(MLL)
matrix XR = (2.3222278923, 2.6066294426, -1.8832365523, -1.3145112385, 1.2728299833, 1.1697045224, 0.3808007836, -0.6188382957)
assert reldif(e(ll), -33.1426170323) < 1e-5
forvalues j = 1/8 {
	assert abs(el(e(b),1,`j') - el(XR,1,`j')) < 1e-2 * (1 + abs(el(XR,1,`j')))
}
forvalues s = 1/5 {
	nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(`s')
	assert reldif(e(ll), scalar(MLL)) < 1e-7
	assert mreldif(e(b), M11) < 1e-3
}
* the seed is set before estimation: the same RNG state afterwards (set seed then fit = seed())
set seed 11
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver
assert mreldif(e(b), M11) == 0
di as text "  PASS: multi engine reproducible with seed(), same optimum for seeds 1-5, equal to goldfish"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: documented e() entries"
di as text "{hline 60}"
nwdynam w40, seed(1)
assert "`e(cmd)'" == "nwdynam"
assert "`e(depvar)'" == "w40"
assert "`e(submodel)'" == "choice"
assert "`e(effects)'" == "inertia recip indeg"
assert "`e(ties)'" == "none"
assert "`e(unidentified)'" == ""
assert "`e(title)'" == "Dynamic Network Actor Model, choice sub-model (conditional logit, MLE)"
assert "`e(version)'" != ""
assert e(N) == 40 & e(nodes) == 4 & e(ntimes) == 40
assert reldif(e(ll), -38.1300248) < 1e-5
assert "`: colnames e(b)'" == "inertia recip indeg"
assert "`: rownames e(V)'" == "inertia recip indeg"
nwdynam w40, submodel(rate) seed(1)
assert "`e(submodel)'" == "rate"
assert "`e(title)'" == "Dynamic Network Actor Model, rate sub-model, no intercept (ordinal partial likelihood, MLE)"
nwdynam w40, submodel(rate) intercept seed(1)
assert "`e(title)'" == "Dynamic Network Actor Model, rate sub-model, WITH intercept (continuous-time competing-risks hazard, MLE)"
* e(effects) lists the e(b) columns, Intercept included (help: e(b) columns named per e(effects))
assert "`e(effects)'" == "Intercept indeg outdeg"
assert "`: colnames e(b)'" == "Intercept indeg outdeg"
* tied times: e(ties) and e(ntimes)
nwclear
xt_w40
replace t = ceil(t / 4)
nwset sender receiver, eventtime(t) directed name(wq)
nwdynam wq, seed(1)
assert "`e(ties)'" == "simultaneous" & e(ntimes) == 10 & e(N) == 40
nwdynam wq, sequential seed(1)
assert "`e(ties)'" == "sequential" & e(ntimes) == 10
* sequential on distinct times: e(ties) stays none
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
nwdynam w40, sequential seed(1)
assert "`e(ties)'" == "none"
assert mreldif(e(b), C1) == 0
di as text "  PASS: e() entries"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: nwdynam version"
di as text "{hline 60}"
nwdynam w40, seed(1)
local ev "`e(version)'"
matrix XB = e(b)
nwdynam version
assert "`r(reference)'" == "goldfish 1.6.12"
assert "`r(license)'" == "GPL-3"
assert "`r(version)'" == "`ev'"
assert `"`r(core)'"' != ""
assert `"`r(plugin)'"' != ""
* the version subcommand leaves the estimation results alone
assert "`e(cmd)'" == "nwdynam" & mreldif(e(b), XB) == 0
di as text "  PASS: nwdynam version"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: postestimation on every result shape"
di as text "{hline 60}"
nwdynam w40, seed(1)
xt_post "choice, native engine"
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(11)
xt_post "choice, multi engine, 8 effects"
nwdynam w40, submodel(rate) seed(1)
xt_post "rate"
nwdynam w40, submodel(rate) intercept seed(1)
xt_post "rate with intercept"
lincom Intercept
assert abs(r(estimate) - (-1.259479036)) < 2e-3
nwclear
xt_w40
nwset sender receiver, eventtime(t) undirected name(wu)
nwdynam wu, submodel(choice_coordination) inertia indeg seed(1)
assert "`e(submodel)'" == "choice_coordination"
assert "`e(title)'" == "Dynamic Network Actor Model, choice_coordination sub-model (multinomial-multinomial, MLE)"
xt_post "choice_coordination"
nwclear
import delimited "dev/dynam_crosscheck_twomode_events.csv", clear varnames(1)
* numeric codes: string ids P../O.. hit CONFIRMED BUG N3 (core: mode split does not match the node order)
encode sender, gen(sn)
encode receiver, gen(rn)
nwset sn rn, twomode eventtime(time) name(tm)
nwdynam tm, seed(1)
assert e(nodes) == 10
xt_post "two-mode choice"
nwdynam tm, submodel(rate) seed(1)
assert "`e(effects)'" == "outdeg"
xt_post "two-mode rate"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block E: a fit leaves data, frame, event list and name space unchanged"
di as text "{hline 60}"
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
nwload w40, xvars
gen double z = _n^2
gen double g = mod(_n, 2)
tempfile snap
quietly save `snap'
local fr "`c(frame)'"
_nwsyntax w40, max(1)
mata: XT_EL0 = *(`netobj'->get_eventlist())
scalar xt_keep_scalar = 17
matrix xt_keep_matrix = (1, 2)
matrix xt_keep_ref = (1, 2)
foreach call in "w40, seed(1)" "w40, inertia recip indeg outdeg trans seed(1)" "w40, same(g) alter(z) seed(1)" ///
	"w40, submodel(rate) seed(1)" "w40, submodel(rate) ego(z) tertius(g) seed(1)" "w40, submodel(rate) intercept seed(1)" ///
	"w40, inertiawindow(5) seed(1)" "w40, weightedinertia seed(1)" {
	nwdynam `call'
	cf _all using `snap'
	assert "`c(frame)'" == "`fr'"
	capture frame __nwdynam_native: describe
	assert _rc != 0
	_nwsyntax w40, max(1)
	mata: st_numscalar("xt_eq", XT_EL0 == *(`netobj'->get_eventlist()))
	assert scalar(xt_eq) == 1
	capture confirm scalar __nwdynam_sequential
	assert _rc != 0
	capture confirm matrix __nwdynam_b
	assert _rc != 0
	capture confirm scalar __nwdynam_ll
	assert _rc != 0
	assert scalar(xt_keep_scalar) == 17
	assert mreldif(xt_keep_matrix, xt_keep_ref) == 0
	di as text "  PASS: nwdynam `call' leaves state unchanged"
}

di "ALL PASS"
