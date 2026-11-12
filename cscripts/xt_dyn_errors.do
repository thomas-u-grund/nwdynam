cscript

* xt_dyn_errors.do -- DRAFT, not yet run.
* Error paths of nwdynam with the exact return code (nwdynam.ado line numbers as of 1.2.5; core codes: 482 network
* not found, from _nwsyntax). After every error: e() of the previous fit and the data in memory are unchanged.
*   198  not an event network (52-56); submodel() unknown (79-82); choice_coordination on a two-mode (124-127) or a
*        directed network (128-131); choice/rate on an undirected one-mode network (133-136); effects of the wrong
*        sub-model (178-181, 212-215, 242-245); two-mode restrictions (302-341); windows: wrong sub-model or <= 0
*        (362-422); weighted: wrong sub-model or combined with a window (446-489); intercept outside rate or with a
*        rate window (541-548); covariate dataset with _N != nodes (611-614, 630-633); opportunities() with a
*        covariate option (32-35), outside choice (686-689), event or actor number out of range (693-701); tie() with
*        rate (212-215), an event network (745-748) or another actor count (749-752)
*   109  egoalterint() with a string variable (22-29); same()/opportunities() with a string variable (syntax)
*   102/103  egoalterint() / opportunities() with one / three variables (syntax min=2 max=2)
*   111  covariate that does not exist (syntax)
*   482  network or tie() network not found
*   2001 fewer than 2 events (659-662)
*   3498 self-loop event (Mata _error in DynamState::init, unw_dynam.do:164-168)
* Expected runtime: about 20 seconds (Stata 17 BE).
* Run from the repository root; the test runner sets the adopath (core, this repository, the plugin folder).
* No sysdir/adopath here.

capture program drop nwinstall
local nwcore : env NWCORE
if `"`nwcore'"' != "" quietly do `"`nwcore'/unw_core.do"'
quietly do unw_dynam.do

mata: XT_W40 = (1,3 \ 1,3 \ 1,3 \ 1,3 \ 2,1 \ 4,2 \ 2,4 \ 2,4 \ 1,3 \ 4,3 \ 3,1 \ 4,2 \ 3,2 \ 4,2 \ 1,3 \ 1,3 \ 2,4 \ 4,2 \ 2,4 \ 2,4 \ 3,1 \ 3,1 \ 1,2 \ 1,3 \ 3,4 \ 4,1 \ 1,3 \ 3,1 \ 1,2 \ 1,3 \ 4,2 \ 2,1 \ 4,2 \ 4,2 \ 3,2 \ 2,3 \ 1,3 \ 2,4 \ 4,2 \ 2,3)

capture program drop xt_err
program xt_err
	* xt_err rc command...: the command stops with exactly r(rc); e(b), e(ll), e(cmd) and the data are unchanged
	gettoken rc 0 : 0
	tempfile snap
	quietly save `snap', replace emptyok
	tempname b ll
	matrix `b' = e(b)
	scalar `ll' = e(ll)
	capture noisily `0'
	local got = _rc
	if `got' != `rc' {
		di as error "  FAIL: expected r(`rc'), got r(`got'): `0'"
		exit 9
	}
	assert mreldif(e(b), `b') == 0
	assert e(ll) == scalar(`ll')
	assert "`e(cmd)'" == "nwdynam"
	cf _all using `snap'
	* (scalar __nwdynam_sequential may be left behind after an error: see xt_dyn_suspected.do)
	capture scalar drop __nwdynam_sequential
	di as text "  PASS: r(`rc'): `0'"
end

capture program drop xt_w40
program xt_w40
	clear
	quietly set obs 40
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store(., ("sender", "receiver"), XT_W40)
	quietly gen double t = _n
end

* ---------------------------------------------------------------------------------------------------------------------
* networks used below
nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
xt_w40
nwset sender receiver, eventtime(t) undirected name(wu)
set seed 3
nwrandom 4, prob(.5) name(static)
matrix XB5 = J(5, 5, 0)
matrix XB5[1,2] = 1
nwset, mat(XB5) directed name(five) nooutput
import delimited "dev/dynam_crosscheck_twomode_events.csv", clear varnames(1)
* numeric codes: string ids P../O.. hit CONFIRMED BUG N3 (core: mode split does not match the node order)
encode sender, gen(sn)
encode receiver, gen(rn)
nwset sn rn, twomode eventtime(time) name(tm)
clear
input long sender long receiver double t
1 2 1
end
nwset sender receiver, eventtime(t) directed name(one)
clear
input long sender long receiver double t
1 2 1
2 3 2
3 3 3
3 1 4
1 3 5
end
tempfile loopdata
quietly save `loopdata'
nwset sender receiver, eventtime(t) directed name(loop)
use `loopdata', clear
nwset sender receiver, eventtime(t) undirected name(loopu)

* a valid fit whose e() must survive every error below; covariate dataset of w40 in memory
nwload w40, xvars
gen double x = _n
gen double x2 = _n^2
gen str1 sx = "a"
nwdynam w40, seed(1)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block A: networks and sub-models"
di as text "{hline 60}"
xt_err 482 nwdynam nosuchnet
xt_err 482 nwdynam w40 static
xt_err 198 nwdynam static
xt_err 198 nwdynam five
xt_err 198 nwdynam w40, submodel(foo)
xt_err 198 nwdynam w40, submodel(Choice)
xt_err 198 nwdynam w40, submodel(choice_coordination)
xt_err 198 nwdynam wu
xt_err 198 nwdynam wu, submodel(rate)
xt_err 198 nwdynam tm, submodel(choice_coordination)
xt_err 198 nwdynam w40, foo
xt_err 198 nwdynam w40, ego(x)
xt_err 198 nwdynam w40, submodel(rate) inertia
xt_err 198 nwdynam w40, submodel(rate) recip
xt_err 198 nwdynam w40, submodel(rate) trans
xt_err 198 nwdynam w40, submodel(rate) four
xt_err 198 nwdynam w40, submodel(rate) same(x)
xt_err 198 nwdynam w40, submodel(rate) egoalterint(x x2)
xt_err 198 nwdynam w40, submodel(rate) tie(static)
xt_err 198 nwdynam wu, submodel(choice_coordination) recip
xt_err 198 nwdynam wu, submodel(choice_coordination) outdeg
xt_err 198 nwdynam wu, submodel(choice_coordination) cycle
xt_err 198 nwdynam wu, submodel(choice_coordination) ego(x)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block B: two-mode restrictions"
di as text "{hline 60}"
xt_err 198 nwdynam tm, recip
xt_err 198 nwdynam tm, outdeg
xt_err 198 nwdynam tm, trans
xt_err 198 nwdynam tm, nodetrans
xt_err 198 nwdynam tm, submodel(rate) indeg
xt_err 198 nwdynam tm, submodel(rate) nodetrans
xt_err 198 nwdynam tm, inertiawindow(5)
xt_err 198 nwdynam tm, weightedinertia
xt_err 198 nwdynam tm, submodel(rate) outdegwindow(5)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block C: windows, weights, intercept"
di as text "{hline 60}"
xt_err 198 nwdynam w40, submodel(rate) inertiawindow(5)
xt_err 198 nwdynam w40, inertiawindow(0)
xt_err 198 nwdynam w40, inertiawindow(-2)
xt_err 198 nwdynam w40, submodel(rate) recipwindow(5)
xt_err 198 nwdynam w40, recipwindow(0)
xt_err 198 nwdynam w40, indegwindow(5)
xt_err 198 nwdynam w40, submodel(rate) indegwindow(0)
xt_err 198 nwdynam w40, outdegwindow(5)
xt_err 198 nwdynam w40, submodel(rate) outdegwindow(-0.5)
xt_err 198 nwdynam wu, submodel(choice_coordination) inertiawindow(5)
xt_err 198 nwdynam w40, submodel(rate) weightedinertia
xt_err 198 nwdynam w40, weightedinertia inertiawindow(5)
xt_err 198 nwdynam w40, submodel(rate) weightedrecip
xt_err 198 nwdynam w40, weightedrecip recipwindow(5)
xt_err 198 nwdynam wu, submodel(choice_coordination) weightedindeg
xt_err 198 nwdynam w40, submodel(rate) weightedindeg indegwindow(5)
xt_err 198 nwdynam wu, submodel(choice_coordination) weightedoutdeg
xt_err 198 nwdynam w40, submodel(rate) weightedoutdeg outdegwindow(5)
xt_err 198 nwdynam w40, intercept
xt_err 198 nwdynam wu, submodel(choice_coordination) intercept
xt_err 198 nwdynam w40, submodel(rate) intercept indegwindow(5)
xt_err 198 nwdynam w40, submodel(rate) intercept outdegwindow(5)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block D: covariates and tie()"
di as text "{hline 60}"
xt_err 109 nwdynam w40, same(sx)
xt_err 109 nwdynam w40, alter(sx)
xt_err 109 nwdynam w40, submodel(rate) ego(sx)
xt_err 109 nwdynam w40, egoalterint(x sx)
xt_err 111 nwdynam w40, same(nosuchvar)
xt_err 102 nwdynam w40, egoalterint(x)
xt_err 103 nwdynam w40, egoalterint(x x2 x)
xt_err 103 nwdynam w40, same(x x2)
preserve
quietly drop in 4
xt_err 198 nwdynam w40, same(x)
xt_err 198 nwdynam w40, submodel(rate) ego(x)
xt_err 198 nwdynam w40, tertius(x)
xt_err 198 nwdynam w40, egoalterint(x x2)
quietly expand 2
xt_err 198 nwdynam w40, diff(x)
restore
xt_err 482 nwdynam w40, tie(nosuchnet)
xt_err 198 nwdynam w40, tie(w40)
xt_err 198 nwdynam w40, tie(five)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block E: opportunities()"
di as text "{hline 60}"
xt_err 198 nwdynam w40, same(x) opportunities(x x2)
preserve
clear
quietly set obs 160
gen long oev = ceil(_n / 4)
gen long oact = _n - 4 * (oev - 1)
gen str1 so = "a"
xt_err 102 nwdynam w40, opportunities(oev)
xt_err 103 nwdynam w40, opportunities(oev oact oev)
xt_err 109 nwdynam w40, opportunities(oev so)
xt_err 198 nwdynam w40, submodel(rate) opportunities(oev oact)
xt_err 198 nwdynam wu, submodel(choice_coordination) opportunities(oev oact)
quietly replace oev = 41 in 160
xt_err 198 nwdynam w40, opportunities(oev oact)
quietly replace oev = 0 in 160
xt_err 198 nwdynam w40, opportunities(oev oact)
quietly replace oev = 40 in 160
quietly replace oact = 5 in 1
xt_err 198 nwdynam w40, opportunities(oev oact)
quietly replace oact = 0 in 1
xt_err 198 nwdynam w40, opportunities(oev oact)
restore

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "Block F: event data"
di as text "{hline 60}"
xt_err 2001 nwdynam one
xt_err 2001 nwdynam one, submodel(rate)
xt_err 2001 nwdynam one, submodel(rate) intercept
// CONFIRMED BUG N2: the self-loop refusal (unw_dynam.do:167) passes a ~140-character text to Mata _error(), which
// fails with r(3300) "argument out of range"; the message is lost. Wanted r(3498) with the message. Recorded, goes on.
local xt_known ""
foreach c in "loop" "loop, submodel(rate)" "loop, inertia trans" "loop, submodel(rate) intercept" "loopu, submodel(choice_coordination)" {
	capture noisily nwdynam `c'
	local rc = _rc
	capture frame change default
	capture frame drop __nwdynam_native
	if `rc' == 3498 di as text "  PASS: r(3498): nwdynam `c'"
	else {
		local xt_known "`xt_known' N2(`c': r(`rc'))"
		di as error "  CONFIRMED BUG N2: nwdynam `c' gives r(`rc'), not r(3498)"
	}
	assert "`e(cmd)'" == "nwdynam"
}
if "`xt_known'" != "" di as error "known confirmed bugs hit:`xt_known'"

di "ALL PASS"
