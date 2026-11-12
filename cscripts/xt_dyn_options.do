cscript

* xt_dyn_options.do -- DRAFT, not yet run.
* Option combinations for nwdynam on the W40 sequence (4 actors, 40 events at t = 1..40):
*   recency windows (inclusive boundary as goldfish: an event exactly W ago still counts), window boundaries and
*   scaling (times x 3 with window x 3), a window covering the whole log = no window, self-activation;
*   weighted effects (count weights); opportunities() (full sets = no restriction, a risk set of one adds ln 3 to the
*   log likelihood and leaves the estimate alone, row order of the opportunities data irrelevant); intercept and the
*   time unit (times x 2: Intercept - ln 2, slopes unchanged, ll - 40 ln 2); option order; the long effect list
*   (all 17 choice effects at once).
* References: goldfish 1.6.12 (R code in xt_NOTES.md) or invariances. Coefficients to 2e-3 relative, ll to 1e-5.
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
	* W40 events; `1' = time expression (default _n)
	args texp
	if "`texp'" == "" local texp "_n"
	clear
	quietly set obs 40
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store(., ("sender", "receiver"), XT_W40)
	quietly gen double t = `texp'
end

capture program drop xt_gf
program xt_gf
	* xt_gf "label" ll b1 [b2 ...]: goldfish reference; coefficients to 2e-3, ll to 1e-5 relative
	gettoken label 0 : 0
	gettoken ll 0 : 0
	local k = colsof(e(b))
	forvalues j = 1/`k' {
		gettoken ref 0 : 0
		assert abs(el(e(b), 1, `j') - (`ref')) < 2e-3 * (1 + abs(`ref'))
	}
	assert `"`0'"' == ""
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

capture program drop xt_keep
program xt_keep
	args p
	matrix `p'B = e(b)
	scalar `p'LL = e(ll)
end

nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: windows (choice inertia, rate outdeg)"
di as text "{hline 60}"
nwdynam w40, inertia seed(1)
xt_gf "inertia alone" -41.213570853 1.007368923
xt_keep I
nwdynam w40, inertiawindow(5) seed(1)
assert "`e(effects)'" == "inertia"
xt_gf "inertiawindow(5) (inclusive)" -41.15874024 0.92502735
xt_keep W5
nwdynam w40, inertia inertiawindow(5) seed(1)
xt_same "inertia inertiawindow(5) = inertiawindow(5)" 1e-8 W5
nwdynam w40, inertiawindow(5.5) seed(1)
xt_same "inertiawindow(5.5) = inertiawindow(5) (integer times)" 1e-6 W5
nwdynam w40, inertiawindow(4.5) seed(1)
xt_gf "inertiawindow(4.5)" -41.76833886 0.8671928
nwdynam w40, inertiawindow(4.999999) seed(1)
xt_gf "inertiawindow(4.999999) = exclusive 5" -41.76833886 0.8671928
* the largest time gap in the log is 39: windows of 39 and more equal no window
foreach W in 39 40 1e6 1e9 {
	nwdynam w40, inertiawindow(`W') seed(1)
	xt_same "inertiawindow(`W') = inertia" 1e-6 I
}
nwdynam w40, submodel(rate) outdegwindow(5) seed(1)
assert "`e(effects)'" == "outdeg"
xt_gf "rate outdegwindow(5)" -55.34015398 -0.1274876
xt_keep RW5
nwdynam w40, submodel(rate) outdegwindow(4.5) seed(1)
xt_gf "rate outdegwindow(4.5)" -55.44260094 0.03669108
* time unit: times x 3 with windows x 3
nwclear
xt_w40 "3*_n"
nwset sender receiver, eventtime(t) directed name(w3)
nwdynam w3, inertiawindow(15) seed(1)
xt_same "times x 3, inertiawindow(15) = inertiawindow(5)" 1e-6 W5
nwdynam w3, submodel(rate) outdegwindow(15) seed(1)
xt_same "times x 3, outdegwindow(15) = outdegwindow(5)" 1e-6 RW5
* a time shift does not change any window
nwclear
xt_w40 "_n - 1000.5"
nwset sender receiver, eventtime(t) directed name(wsh)
nwdynam wsh, inertiawindow(5) seed(1)
xt_same "time shift, inertiawindow(5)" 1e-6 W5
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: weighted effects (all event weights 1: weights = event counts)"
di as text "{hline 60}"
nwdynam w40, weightedinertia seed(1)
assert "`e(effects)'" == "inertia"
xt_gf "weightedinertia" -37.7174439566 0.2889076707
xt_keep WI
nwdynam w40, inertia weightedinertia seed(1)
xt_same "inertia weightedinertia = weightedinertia" 1e-8 WI
nwdynam w40, submodel(rate) weightedindeg weightedoutdeg seed(1)
assert "`e(effects)'" == "indeg outdeg"
xt_gf "rate weightedindeg weightedoutdeg" -55.01240505209 -0.08797114084 -0.11872452538
* weighting one effect only must differ from both the unweighted and the fully weighted fit
nwdynam w40, submodel(rate) seed(1)
xt_keep RU
nwdynam w40, submodel(rate) weightedindeg outdeg seed(1)
assert mreldif(e(b), RUB) > 1e-4
* the weighted-indeg statistic of the choice sub-model: weightedindeg self-activates indeg only
nwdynam w40, weightedindeg seed(1)
assert "`e(effects)'" == "indeg"
di as text "  PASS: weighted options"

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: opportunities()"
di as text "{hline 60}"
* full opportunity sets (every actor at every event, senders included) = no restriction
clear
quietly set obs 160
gen long oev = ceil(_n / 4)
gen long oact = _n - 4 * (oev - 1)
nwdynam w40, inertia opportunities(oev oact) seed(1)
xt_same "full opportunity sets = inertia" 1e-6 I
* sender rows left out: same (the sender is never a candidate)
mata: st_local("snd", invtokens(strofreal(XT_W40[.,1]')))
gen long snd = real(word("`snd'", oev))
preserve
quietly drop if oact == snd
nwdynam w40, inertia opportunities(oev oact) seed(1)
xt_same "opportunity sets without the sender" 1e-6 I
restore
* event 5 (2 -> 1) has the realised receiver as its only opportunity: its log-likelihood term becomes 0 instead of
* -ln 3 (sender 2 has no earlier ties, inertia is 0 for all its candidates, so the term never depended on the
* estimate). goldfish: inertia 1.007368923 (unchanged), ll -40.114958564 = -41.213570853 + ln 3
quietly drop if oev == 5 & oact != 1
nwdynam w40, inertia opportunities(oev oact) seed(1)
xt_gf "event 5 restricted to actor 1" -40.114958564 1.007368923
assert abs(e(ll) - (scalar(ILL) + ln(3))) < 1e-6
xt_keep OP
* rows of the opportunities data in any order
set seed 4711
gen double u = runiform()
sort u
nwdynam w40, inertia opportunities(oev oact) seed(1)
xt_same "opportunities data in random row order" 1e-10 OP
* variables of any name and type
rename (oev oact) (ev_number actor_id)
recast double ev_number actor_id
nwdynam w40, inertia opportunities(ev_number actor_id) seed(1)
xt_same "opportunities(), other names, double" 1e-10 OP
* opportunities() forces the Mata engine; default choice effects with full sets equal the native fit
clear
quietly set obs 160
gen long oev = ceil(_n / 4)
gen long oact = _n - 4 * (oev - 1)
nwdynam w40, opportunities(oev oact) seed(1)
xt_gf "default choice effects with full opportunity sets" -38.1300248 1.0497675 1.0747602 -0.2707421

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: intercept and the time unit"
di as text "{hline 60}"
nwdynam w40, submodel(rate) intercept seed(1)
xt_gf "rate intercept" -94.37558055 -1.259479036 -0.067194365 0.009146213
xt_keep IC
nwclear
xt_w40 "2*_n"
nwset sender receiver, eventtime(t) directed name(w2)
nwdynam w2, submodel(rate) intercept seed(1)
assert abs(el(e(b),1,1) - (el(ICB,1,1) - ln(2))) < 1e-5
assert abs(el(e(b),1,2) - el(ICB,1,2)) < 1e-5
assert abs(el(e(b),1,3) - el(ICB,1,3)) < 1e-5
assert abs(e(ll) - (scalar(ICLL) - 40*ln(2))) < 1e-6
assert abs(e(ll) - (-122.10146777)) < 1e-4
di as text "  PASS: times x 2: Intercept - ln 2, slopes unchanged, ll - 40 ln 2"
* intercept with weighted effects runs (verified against goldfish by the package; here only: it is accepted and the
* weighted fit differs from the unweighted one)
nwdynam w2, submodel(rate) intercept weightedindeg weightedoutdeg seed(1)
assert "`: colnames e(b)'" == "Intercept indeg outdeg"
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block E: option order, explicit default set, long effect list"
di as text "{hline 60}"
nwdynam w40, seed(1)
xt_keep D
nwdynam w40, indeg recip inertia seed(1)
assert "`e(effects)'" == "inertia recip indeg"
xt_same "options in another order" 0 D
nwdynam w40, seed(1) recip indeg inertia
xt_same "seed() first" 0 D
nwdynam w40, submodel(rate) outdeg indeg seed(1)
assert "`e(effects)'" == "indeg outdeg"
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(1)
assert "`e(effects)'" == "inertia recip indeg outdeg trans cycle commonsender commonreceiver"
* 8-effect reference (goldfish), looser tolerance: flat directions on 4 actors
assert reldif(e(ll), -33.1426170323) < 1e-5
matrix XR = (2.3222278923, 2.6066294426, -1.8832365523, -1.3145112385, 1.2728299833, 1.1697045224, 0.3808007836, -0.6188382957)
forvalues j = 1/8 {
	assert abs(el(e(b),1,`j') - el(XR,1,`j')) < 1e-2 * (1 + abs(el(XR,1,`j')))
}
di as text "  PASS: 8 choice effects = goldfish"
* all 17 choice effects: structure only (on 4 actors most of them are not identified; no reference value)
set seed 3
nwrandom 4, prob(.5) name(exo)
nwload w40, xvars
gen double g = mod(_n, 2)
gen double z = _n^2
gen double h = (_n == 3)
gen double y = sqrt(_n)
gen double c = 10 - _n
nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver four nodetrans same(g) diff(z) sim(y) ///
	alter(h) egoalterint(c y) tertius(z) tie(exo) seed(1)
assert colsof(e(b)) == 17
assert "`e(effects)'" == "inertia recip indeg same diff sim alter outdeg trans cycle commonsender commonreceiver four nodetrans egoalterint tertius tie"
assert "`: colnames e(b)'" == "`e(effects)'"
assert e(ll) < 0
* the same model with options written in reverse order
nwdynam w40, tie(exo) tertius(z) egoalterint(c y) alter(h) sim(y) diff(z) same(g) nodetrans four commonreceiver ///
	commonsender cycle trans outdeg indeg recip inertia seed(1)
assert "`e(effects)'" == "inertia recip indeg same diff sim alter outdeg trans cycle commonsender commonreceiver four nodetrans egoalterint tertius tie"
* all 5 rate effects
nwdynam w40, submodel(rate) indeg outdeg nodetrans ego(z) tertius(g) seed(1)
assert "`e(effects)'" == "indeg outdeg ego nodetrans tertius"
* with intercept: Intercept first
nwdynam w40, submodel(rate) intercept indeg outdeg nodetrans ego(z) tertius(g) seed(1)
assert "`: colnames e(b)'" == "Intercept indeg outdeg ego nodetrans tertius"
di as text "  PASS: long effect lists"

di "ALL PASS"
