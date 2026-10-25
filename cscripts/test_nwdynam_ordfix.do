cscript

* Regression test: events with the same time (ties). The events were sorted by time with order() on the time column alone,
* which breaks ties arbitrarily (it depends on the state of the random-number generator), so the same data gave a different
* order of the tied events, a different tie history and a different likelihood from run to run. Tied events now keep the
* order in which the network lists them. Check: (a) the same fit under different seeds and after random draws gives the
* identical log likelihood; (b) it equals the fit of the same events with strictly increasing times in that order.

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

nwclear
nwwebuse socialevolution, nwclear
generate tt = ceil(t/20)
generate long row = _n
sort tt row
duplicates report tt
generate double tstrict = tt + row/1e6
tempfile d
save `d'

* strictly ordered reference (same event order, no ties)
nwset sender receiver, eventtime(tstrict) directed name(strict)
quietly nwdynam strict
scalar llref = e(ll)
matrix bref = e(b)

use `d', clear
nwset sender receiver, eventtime(tt) directed name(tied)
forvalues s = 1/6 {
	set seed `s'
	forvalues j = 1/`s' {
		quietly mata: junk = runiform(1000, 1)
	}
	quietly nwdynam tied
	assert reldif(e(ll), llref) < 1e-8
	assert mreldif(e(b), bref) < 1e-5
	di as text "  PASS: seed `s': ll " %12.5f e(ll) ", same as the strictly ordered events"
}
di as result "ALL TESTS PASSED"
