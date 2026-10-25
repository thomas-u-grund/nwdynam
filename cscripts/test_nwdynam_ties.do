cscript

* Events with an identical time (ties) are SIMULTANEOUS, as in goldfish: the statistics of every event of a time point are
* computed from the state before that time point; the state is updated after the whole group. The numbers below come from
* goldfish 1.6.12 (dev/dynam_ties_crosscheck.R) on the same events. Also checked: the order of the tied rows does not
* matter, option sequential keeps the old behaviour, the native and the Mata evaluators agree under ties, data without
* ties are unchanged (the unchanged-estimates part is checked in test_nwdynam_ado/native/fixes1).

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
global nwwebpath `"`nwcore'/data"'

* compare e(b) and e(ll) with goldfish: coefficients to 2e-3 (the two optimizers stop at slightly different points), ll to 1e-5 relative
capture program drop chk
program chk
	args label ll
	local k = colsof(e(b))
	forvalues j = 1/`k' {
		local ref ``=`j'+2''
		assert abs(el(e(b), 1, `j') - `ref') < 2e-3
	}
	assert reldif(e(ll), `ll') < 1e-5
	di as text "  PASS: `label': ll " %12.5f e(ll) " (goldfish " %12.5f `ll' ")"
end

nwclear
nwwebuse socialevolution, nwclear
generate tt = ceil(t/20)
generate th = ceil(t*2)
generate long row = _n
tempfile d
save `d'

* ---- A: 4 time points, 435 tied events
use `d', clear
nwset sender receiver, eventtime(tt) directed name(a)
nwdynam a, inertia recip indeg
assert e(ntimes) == 4
assert "`e(ties)'" == "simultaneous"
chk "A choice inertia recip indeg" -1595.32649 3.5380732 1.3654605 -.3843064
nwdynam a, submodel(rate) indeg outdeg
chk "A rate indeg outdeg" -1735.62725 .3342426 .1606922

* ---- B: 47 time points
nwclear
use `d', clear
nwset sender receiver, eventtime(th) directed name(b)
nwdynam b, inertia recip indeg
chk "B choice (native)" -825.301403 4.584526 1.339734 -.358509
scalar llnative = e(ll)
nwdynam b, inertia recip indeg inertiawindow(1e9)
assert reldif(e(ll), llnative) < 1e-9
di as text "  PASS: B choice, Mata evaluator = native evaluator under ties"
nwdynam b, inertia recip indeg outdeg trans
chk "B choice + outdeg trans" -808.384668 4.4871244 1.889513 -.3468104 -.3027759 .5348534
nwdynam b, inertia recip indeg inertiawindow(5)
chk "B choice inertiawindow(5)" -908.015839 4.0422002 2.17348962 -.07505825
nwdynam b, submodel(rate) indeg outdeg
chk "B rate indeg outdeg" -1674.11949 .2603611 .239943
nwdynam b, submodel(rate) intercept indeg outdeg
chk "B rate intercept" -1271.82665 -2.9698256 .5391501 .2448834

* ---- order of the tied rows does not matter
nwclear
use `d', clear
set seed 11
generate double u = runiform()
sort th u
nwset sender receiver, eventtime(th) directed name(b2)
nwdynam b2, inertia recip indeg
chk "B rows of ties shuffled" -825.301403 4.584526 1.339734 -.358509
nwdynam b2, submodel(rate) intercept indeg outdeg
chk "B rate intercept, shuffled" -1271.82665 -2.9698256 .5391501 .2448834

* ---- option sequential: the old behaviour (events one after another, in input order)
nwclear
use `d', clear
sort tt row
nwset sender receiver, eventtime(tt) directed name(c)
nwdynam c, inertia recip indeg sequential
assert "`e(ties)'" == "sequential"
assert reldif(e(ll), -652.8339) < 1e-6
assert abs(el(e(b), 1, 1) - 4.99798) < 1e-4
di as text "  PASS: sequential reproduces the old fit (ll -652.8339, inertia 5.00)"
nwdynam c, inertia recip indeg
assert e(ll) < -1500

* ---- toy data with tied times: undirected coordination and directed models
nwclear
import delimited "dev/dynam_ties_toy_events.csv", clear varnames(1)
gen s = real(substr(sender, 2, .))
gen r = real(substr(receiver, 2, .))
nwset s r, eventtime(time) undirected name(tc)
nwdynam tc, submodel(choice_coordination) inertia
chk "toy coordination" -37.4340124 .3389611
nwclear
import delimited "dev/dynam_ties_toy_events.csv", clear varnames(1)
gen s = real(substr(sender, 2, .))
gen r = real(substr(receiver, 2, .))
nwset s r, eventtime(time) directed name(td)
nwdynam td, inertia recip indeg
chk "toy choice" -20.5896781 1.4920436 -.1801964 -.8036452
nwdynam td, submodel(rate) indeg outdeg
chk "toy rate" -24.1014731 .2938009 .442106
nwdynam td, submodel(rate) intercept indeg outdeg
chk "toy rate intercept" -26.8764901 -.9173658 -.2049949 .1516609

* ---- data without ties: e(ties) is none
nwclear
use `d', clear
nwset sender receiver, eventtime(t) directed name(e)
nwdynam e, inertia recip indeg
assert "`e(ties)'" == "none"
assert e(ntimes) == e(N)
di as result "ALL TESTS PASSED"
