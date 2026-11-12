cscript

* Clear errors: opportunities() with a covariate option, tie()/opportunities() on a two-mode network,
* egoalterint() with a string variable. Run from the repository root with NWCORE set.
capture program drop nwinstall
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`nwcore'/lib"'
adopath ++ `"`c(pwd)'"'
adopath ++ `"`c(pwd)'/lib"'
capture adopath ++ "lib/plugins/macos"
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do

* --- one-mode toy event network
nwclear
clear
input long sender long receiver long eventtime
1 2 1
2 3 2
3 4 3
4 5 4
5 1 5
1 3 6
2 4 7
3 5 8
4 1 9
5 2 10
end
nwset sender receiver, eventtime(eventtime) name(oppsmall)

* opportunities() with a covariate option is refused before any data are read
clear
input long oppevent long oppactor
1 2
1 3
2 4
2 5
end
gen g = 1
gen h = 2
foreach o in "same(g)" "diff(g)" "sim(g)" "alter(g)" "tertius(g)" "egoalterint(g h)" {
	capture noisily nwdynam oppsmall, inertia opportunities(oppevent oppactor) `o'
	assert _rc == 198
}
capture noisily nwdynam oppsmall, submodel(rate) opportunities(oppevent oppactor) ego(g)
assert _rc == 198
* (opportunities() alone is tested with a full opportunity list in test_nwdynam_ado.do)
di "PASS opportunities with covariates"

* egoalterint() with a string variable: a message that names the option
clear
set obs 5
gen str1 s = "a"
gen a = mod(_n, 2)
gen b = _n
capture noisily nwdynam oppsmall, egoalterint(s a)
assert _rc == 109
capture noisily nwdynam oppsmall, egoalterint(a s)
assert _rc == 109
nwdynam oppsmall, egoalterint(a b)
assert _rc == 0
assert "`e(effects)'" == "egoalterint"
di "PASS egoalterint string"

* --- two-mode event network: tie() and opportunities() are refused
nwclear
clear
import delimited "dev/dynam_crosscheck_twomode_events.csv", clear varnames(1)
encode sender, gen(sendernum)
encode receiver, gen(receivernum)
gen eventtime = time
nwset sendernum receivernum, twomode eventtime(eventtime) name(afftest)
nwdynam afftest
assert _rc == 0
nwsummarize afftest
local n = r(nodes)
mata: __fx_T = J(`n', `n', 0)
nwset, mat(__fx_T) name(covnet) directed nooutput
capture noisily nwdynam afftest, tie(covnet)
assert _rc == 198
clear
input long oppevent long oppactor
1 7
2 8
end
capture noisily nwdynam afftest, opportunities(oppevent oppactor)
assert _rc == 198
mata: mata drop __fx_T
di "PASS two-mode tie opportunities"

di "ALL PASS"
