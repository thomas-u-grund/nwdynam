cscript

* Combinations of nwdynam options against goldfish 1.6.12 (R 4.6): (a) windowed rate effects with intercept, (b) tertiuswindow(),
* (c) inertia/indeg windows and weights in choice_coordination, (d) attribute and closure effects in two-mode networks,
* (e) opportunities() together with attribute effects (oppframe()), (f) change() in two-mode networks and with opportunities().
* The logs and goldfish's estimates, standard errors and log likelihoods are the files of dev/win (made by
* dev/dynam_win_crosscheck.R; all event times are multiples of 0.25, see GOTCHA.md). Tolerance: estimates, standard errors and ll to
* 1e-6 (relative); choice_coordination: estimates 5e-4, standard errors 10% (flat, non-concave likelihood; numerical derivative).
* Then fixed reference values and the error paths of what stays refused.
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
global WINDIR `"`c(pwd)'/dev/win"'
global WINSEC ""
global WINONLY ""
quietly do dev/dynam_win_crosscheck.do
runall

* ---- fixed reference values (goldfish 1.6.12)
quietly runcase a01 rate "indegwindow(6) intercept" one "" 0 0
assert reldif(el(e(b), 1, 1), -2.12372930395133) < 1e-6 & reldif(el(e(b), 1, 2), -0.771118952375098) < 1e-6 & reldif(e(ll), -166.211011578508) < 1e-6
quietly runcase b02 rate "tertius(x) tertiuswindow(6)" one "" 0 0
assert reldif(el(e(b), 1, 1), -0.535670221886186) < 1e-6 & reldif(e(ll), -115.283954858049) < 1e-6
quietly runcase d01 choice "alter(x)" two "" 0 0
assert reldif(el(e(b), 1, 1), -0.0922498017025478) < 1e-6 & reldif(e(ll), -96.2470789948496) < 1e-6
quietly runcase d02 rate "ego(x)" two "" 0 0
assert reldif(el(e(b), 1, 1), -0.0985551504446989) < 1e-6 & reldif(e(ll), -107.235180270851) < 1e-6
quietly runcase e01 choice "same(g)" one "" 1 0
assert reldif(el(e(b), 1, 1), -0.24061711729447) < 1e-6 & reldif(e(ll), -92.6241325197657) < 1e-6
quietly runcase f61 choice "inertia indeg alter(x)" one "x" 1 0
assert reldif(el(e(b), 1, 3), 0.0133680824077382) < 1e-6 & reldif(e(ll), -135.690384258468) < 1e-6
di as result "PASS: fixed reference values"

* ---- error paths: what stays refused (all give error 198; the missing frame gives 111)
foreach spec in ///
	"b02|rate|tertiuswindow(5)|one|-|0" ///
	"b02|choice|tertius(x) tertiuswindow(0)|one|-|0" ///
	"d01|choice|trans|two|-|0" ///
	"d01|choice|cycle|two|-|0" ///
	"d01|choice|commonsender|two|-|0" ///
	"d01|choice|commonreceiver|two|-|0" ///
	"d01|choice|nodetrans|two|-|0" ///
	"d01|choice|same(g)|two|-|0" ///
	"d01|choice|diff(x)|two|-|0" ///
	"d01|choice|sim(x)|two|-|0" ///
	"d01|choice|recip|two|-|0" ///
	"d01|choice|outdeg|two|-|0" ///
	"d01|choice|inertiawindow(3)|two|-|0" ///
	"d01|choice|weightedinertia|two|-|0" ///
	"d02|rate|indeg|two|-|0" ///
	"d02|rate|tertius(x)|two|-|0" ///
	"d02|rate|nodetrans|two|-|0" ///
	"c01|choice_coordination|inertiawindow(3) tertiuswindow(2)|und|-|0" {
	tokenize `"`spec'"', parse("|")
	local chg "`9'"
	if "`chg'" == "-" local chg ""
	capture runcase `1' `3' "`5'" `7' "`chg'" `11' 0
	if _rc != 198 {
		di as error "`spec': expected r(198), got r(" _rc ")"
		exit 9
	}
}
* a change of the tertius variable together with tertiuswindow() and the intercept in a rate model
capture runcase b34 rate "tertius(x) tertiuswindow(4.5) intercept" one "x" 0 0
assert _rc == 198
* opportunities() with a covariate option or change() needs oppframe(); oppframe() needs opportunities() and an existing frame
global WINNOFRAME 1
capture runcase e01 choice "same(g)" one "" 1 0
assert _rc == 198
global WINNOFRAME 2
capture runcase e01 choice "same(g)" one "" 1 0
assert _rc == 198
global WINNOFRAME 3
capture runcase e01 choice "same(g)" one "" 1 0
assert _rc == 111
global WINNOFRAME ""
di as result "PASS: error paths (what goldfish refuses or does not give stays refused)"
