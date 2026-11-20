cscript

* Plot options of nwdynam on the common convention of the nw packages (graph helper _nw_gfinish of the core, 3.29.0 or later):
* estat outliers / changepoints (released plots: default graph names and PNG files unchanged; saving(), export() and twoway options added),
* estat choiceplot / rateplot (default names choiceplot and rateplot, black default titles, title("") suppresses, saving(), export(), noplot).
* Checks: the graph exists, files are written, name()/saving()/export()/graph options without a plot give r(198) (name() alone stays accepted for
* outliers/changepoints), bad values give r(198), the data in memory is untouched, title("") runs, the old-core message (r(199)), and the default
* pictures of the released plots are byte-identical to the released package (reference directory in NWDYNAM_PLOT_REF, made with
* cscripts/nwdynam_plots_defaults.do and the released package; skipped when not set).
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
set graphics off
local pwd `"`c(pwd)'"'

tempfile tmpd
local tmp `"`tmpd'_nwdplots"'
mkdir `"`tmp'"'

* ---- 0. default pictures of the released plots: identical PNG files
local ref : env NWDYNAM_PLOT_REF
if `"`ref'"' == "" {
	di as text "  SKIPPED: default pictures vs released package (set NWDYNAM_PLOT_REF)"
}
else {
	global NWP_OUT `"`tmp'/defaults"'
	quietly do `"`pwd'/cscripts/nwdynam_plots_defaults.do"'
	foreach f in outliers changepoints outliers_named {
		quietly checksum `"`ref'/`f'.png"'
		local c0 = r(checksum)
		local s0 = r(filelen)
		quietly checksum `"$NWP_OUT/`f'.png"'
		if `c0' != r(checksum) | `s0' != r(filelen) {
			di as error "default picture `f'.png differs from the released package"
			exit 9
		}
	}
	di as result "PASS default pictures of estat outliers / changepoints are byte-identical to the released package"
}

* ---- toy data and fits
nwclear
set seed 20261120
clear
set obs 60
generate long sender = ceil(6 * runiform())
generate long receiver = ceil(5 * runiform())
replace receiver = receiver + (receiver >= sender)
generate double eventtime = _n + runiform() / 4
sort eventtime
nwset sender receiver, eventtime(eventtime) name(toy) directed
generate double canary = _n^2
datasignature set, reset
local sig = r(datasignature)
nwdynam toy, submodel(choice) inertia

* ---- 1. estat outliers / changepoints
graph drop _all
estat outliers, method(iqr) parameter(1) plot
assert r(n_outliers) >= 0 & r(N) == e(N)
graph describe Graph
graph drop _all
estat outliers, method(iqr) parameter(1) plot name(go1) saving("`tmp'/ol1") export("`tmp'/ol1.png", width(500)) title("Own title") ylabel(, angle(0))
graph describe go1
confirm file `"`tmp'/ol1.gph"'
confirm file `"`tmp'/ol1.png"'
graph drop _all
estat outliers, method(iqr) parameter(1) plot name(go2, replace) saving("`tmp'/ol2.gph", replace)
graph describe go2
confirm file `"`tmp'/ol2.gph"'
* without a name: the unnamed graph is saved and exported
estat outliers, method(iqr) parameter(1) plot saving("`tmp'/ol3") export("`tmp'/ol3.pdf")
confirm file `"`tmp'/ol3.gph"'
confirm file `"`tmp'/ol3.pdf"'
estat changepoints, window(3) plot name(gc1) saving("`tmp'/cp1") export("`tmp'/cp1.png") title("")
graph describe gc1
confirm file `"`tmp'/cp1.gph"'
confirm file `"`tmp'/cp1.png"'
estat changepoints, window(3) plot
graph describe Graph
* options without plot
foreach c in outliers changepoints {
	capture estat `c', saving("`tmp'/x")
	assert _rc == 198
	capture estat `c', export("`tmp'/x.png")
	assert _rc == 198
	capture estat `c', title("T")
	assert _rc == 198
	capture estat `c', ylabel(0(1)2)
	assert _rc == 198
	capture estat `c', plot saving("`tmp'/x", append)
	assert _rc == 198
	capture estat `c', plot name(bad name)
	assert _rc == 198
	capture estat `c', plot export("`tmp'/nodir/x.png")
	assert _rc == 603
}
assert !fileexists(`"`tmp'/x.gph"') & !fileexists(`"`tmp'/x.png"')
* name() alone without plot stays accepted (it was before)
estat outliers, method(iqr) parameter(1) name(gz)
estat changepoints, window(3) name(gz)
* a bad value is still reported first
capture estat outliers, method(xyz) plot
assert _rc == 198
di as result "PASS estat outliers / changepoints: names, saving, export, pass-through, errors"

* ---- 2. estat choiceplot
graph drop _all
estat choiceplot, topk(1 3 5)
assert r(N) == e(N)
assert r(mean_riskset) == 5
graph describe choiceplot
estat choiceplot, name(gch) saving("`tmp'/ch") export("`tmp'/ch.png", width(500)) xlabel(1(1)5)
graph describe gch
confirm file `"`tmp'/ch.gph"'
confirm file `"`tmp'/ch.png"'
graph drop _all
estat choiceplot, title("") name(gch2, replace) saving("`tmp'/ch2.gph", replace)
graph describe gch2
confirm file `"`tmp'/ch2.gph"'
estat choiceplot, title("My title", color(red))
graph describe choiceplot
* noplot with plot options, bad values
capture estat choiceplot, noplot name(gx)
assert _rc == 198
capture estat choiceplot, noplot saving("`tmp'/x")
assert _rc == 198
capture estat choiceplot, noplot export("`tmp'/x.png")
assert _rc == 198
capture estat choiceplot, noplot title("T")
assert _rc == 198
capture estat choiceplot, name(bad name)
assert _rc == 198
capture estat choiceplot, export("`tmp'/nodir/x.png")
assert _rc == 603
capture estat choiceplot, maxrank(0)
assert _rc == 198
capture estat choiceplot, topk(0)
assert _rc != 0
assert !fileexists(`"`tmp'/x.gph"') & !fileexists(`"`tmp'/x.png"')
estat choiceplot, noplot frame(chf)
assert r(N) == e(N)
di as result "PASS estat choiceplot: names, saving, export, pass-through, errors"

* ---- 3. estat rateplot
nwdynam toy, submodel(rate) indeg
graph drop _all
estat rateplot
assert r(N_actors) == 6
graph describe rateplot
estat rateplot, over(time) periods(5) topactors(3) name(grt) saving("`tmp'/rt") export("`tmp'/rt.png") ylabel(, angle(0))
graph describe grt
confirm file `"`tmp'/rt.gph"'
confirm file `"`tmp'/rt.png"'
graph drop _all
estat rateplot, over(time) periods(5) title("")
graph describe rateplot
estat rateplot, name(gra) saving("`tmp'/ra") export("`tmp'/ra.png") title("")
graph describe gra
confirm file `"`tmp'/ra.gph"'
confirm file `"`tmp'/ra.png"'
capture estat rateplot, noplot name(gx)
assert _rc == 198
capture estat rateplot, noplot saving("`tmp'/x")
assert _rc == 198
capture estat rateplot, noplot ylabel(0(5)20)
assert _rc == 198
capture estat rateplot, name(bad name)
assert _rc == 198
capture estat rateplot, export("`tmp'/nodir/x.png")
assert _rc == 603
capture estat rateplot, over(foo)
assert _rc == 198
estat rateplot, noplot
assert r(N_actors) == 6
assert !fileexists(`"`tmp'/x.gph"')
di as result "PASS estat rateplot: names, saving, export, pass-through, errors"

* ---- 4. the data in memory is untouched, return values survive the graph code
datasignature
assert r(datasignature) == "`sig'"
estat outliers, method(iqr) parameter(1) plot saving("`tmp'/ol9")
nwdynam toy, submodel(choice) inertia
estat outliers, method(iqr) parameter(1) plot saving("`tmp'/ol9")
assert r(n_outliers) < . & r(N) == e(N)
estat changepoints, window(3) plot export("`tmp'/cp9.png")
assert r(n_cpts) < . & r(N) == e(N)
estat choiceplot, saving("`tmp'/ch9")
assert r(N) == e(N) & rowsof(r(curve)) == 20
datasignature
assert r(datasignature) == "`sig'"
di as result "PASS data untouched, returns kept"

* ---- 5. core too old: a clear message and r(199) at plot time only
capture program drop _nw_gfinish
program define _nw_gfinish
	exit 199
end
capture estat choiceplot
assert _rc == 199
capture estat outliers, method(iqr) parameter(1) plot
assert _rc == 199
estat outliers, method(iqr) parameter(1)
estat choiceplot, noplot
program drop _nw_gfinish
estat choiceplot
graph describe choiceplot
di as result "PASS an old core stops plots with r(199); tables unaffected"

di as result "ALL TESTS PASSED (nwdynam plots)"
