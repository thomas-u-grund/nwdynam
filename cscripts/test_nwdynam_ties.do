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
nwwebuse manufacturing, nwclear
keep if t < 5
generate tt = ceil(t*2)
generate th = ceil(t*48)
assert t == t[_n-1] | t > t[_n-1] if _n > 1
generate long tstrict = _n
generate long row = _n
tempfile d
save `d'

* ---- A: 4 time points, 631 tied events
use `d', clear
nwset sender receiver, eventtime(tt) directed name(a)
nwdynam a, inertia recip indeg
assert e(ntimes) == 4
assert "`e(ties)'" == "simultaneous"
chk "A choice inertia recip indeg" -2576.748 1.457867 2.326090 .05560886
nwdynam a, submodel(rate) indeg outdeg
chk "A rate indeg outdeg" -2884.893 .01549955 .1719255

* ---- B: 43 time points
nwclear
use `d', clear
nwset sender receiver, eventtime(th) directed name(b)
nwdynam b, inertia recip indeg
chk "B choice (native)" -2403.316 1.327586 2.486849 .08555338
scalar llnative = e(ll)
nwdynam b, inertia recip indeg inertiawindow(1e9)
assert reldif(e(ll), llnative) < 1e-9
di as text "  PASS: B choice, Mata evaluator = native evaluator under ties"
nwdynam b, inertia recip indeg outdeg trans
chk "B choice + outdeg trans" -2397.745 1.308292 2.484381 .07508011 .01725865 .3528997
nwdynam b, inertia recip indeg inertiawindow(5)
chk "B choice inertiawindow(5)" -2416.93384 1.1708973 2.9183053 .1022024
nwdynam b, submodel(rate) indeg outdeg
chk "B rate indeg outdeg" -2869.869 .02597912 .1678525
nwdynam b, submodel(rate) intercept indeg outdeg
chk "B rate intercept" -2274.449 -2.706195 -.0003323426 .03730036

* ---- order of the tied rows does not matter
nwclear
use `d', clear
set seed 11
generate double u = runiform()
sort th u
nwset sender receiver, eventtime(th) directed name(b2)
nwdynam b2, inertia recip indeg
chk "B rows of ties shuffled" -2403.316 1.327586 2.486849 .08555338
nwdynam b2, submodel(rate) intercept indeg outdeg
chk "B rate intercept, shuffled" -2274.449 -2.706195 -.0003323426 .03730036

* ---- option sequential: the old behaviour (events one after another, in input order)
nwclear
use `d', clear
sort tt row
nwset sender receiver, eventtime(tt) directed name(c)
nwdynam c, inertia recip indeg sequential
assert "`e(ties)'" == "sequential"
assert reldif(e(ll), -2184.5636) < 1e-6
assert abs(el(e(b), 1, 1) - 1.70191) < 1e-4
di as text "  PASS: sequential reproduces the old fit (ll -2184.5636, inertia 1.70)"
nwdynam c, inertia recip indeg
assert e(ll) < -2500

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
nwset sender receiver, eventtime(tstrict) directed name(e)
nwdynam e, inertia recip indeg
assert "`e(ties)'" == "none"
assert e(ntimes) == e(N)
di as result "ALL TESTS PASSED"
