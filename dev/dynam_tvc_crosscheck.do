* nwdynam time-varying covariates against goldfish: reads the files of dev/dynam_tvc_crosscheck.R (environment variable TVDIR)
* and compares estimates, standard errors and log likelihood (tolerance TVTOL, default 1e-6).
* The caller sets the paths (adopath of the core and of this package, a scratch PLUS directory) and the global TVDIR.
local tol : env TVTOL
if "`tol'" == "" local tol 1e-6

import delimited using `"$TVDIR/results.csv"', clear varnames(1)
tempfile res
save `res'
capture frame drop gf
frame create gf
frame gf: use `res', clear

capture program drop runcase
program runcase
	args id sub opts directed inc
	tempfile nodes
	import delimited using `"$TVDIR/nodes_`id'.csv"', clear varnames(1) stringcols(1)
	save `nodes'
	nwclear
	import delimited using `"$TVDIR/events_`id'.csv"', clear varnames(1) stringcols(2 3)
	if "`directed'" == "0" nwset sender receiver, eventtime(time) undirected name(net)
	else nwset sender receiver, eventtime(time) name(net)
	local chgopt ""
	foreach v in x y g {
		capture confirm file `"$TVDIR/changes_`v'_`id'.csv"'
		if !_rc {
			preserve
			import delimited using `"$TVDIR/changes_`v'_`id'.csv"', clear varnames(1) stringcols(2)
			generate node2 = node
			nwset node node2 value, eventtime(time) selfloop name(chg_`v')
			restore
			local chgopt "`chgopt' `v' chg_`v'"
		}
	}
	* node dataset: one row per actor of net, matched by label
	rename _nwnode label
	merge 1:1 label using `nodes', keep(match master) nogen
	rename label _nwnode
	nwdynam net, submodel(`sub') `opts' change(`chgopt') `=cond("`inc'"=="1","increment","")'
end

capture program drop chkcase
program chkcase
	args id tolb tolse
	local tol : env TVTOL
	if "`tol'" == "" local tol 1e-6
	if "`tolb'" == "" local tolb `tol'
	if "`tolse'" == "" local tolse `tolb'
	frame gf: generate long __rn = _n
	frame gf: quietly summarize __rn if case == "`id'"
	frame gf: local r0 = r(min)
	local k = colsof(e(b))
	local worst 0
	forvalues j = 1/`k' {
		frame gf: local gb = estimate[`=`r0' + `j' - 1']
		frame gf: local gs = se[`=`r0' + `j' - 1']
		local d1 = abs(el(e(b), 1, `j') - `gb') / (1 + abs(`gb'))
		local d2 = abs(sqrt(el(e(V), `j', `j')) - `gs') / (1 + abs(`gs'))
		if `d1' / `tolb' > `worst' local worst = `d1' / `tolb'
		if `d2' / `tolse' > `worst' local worst = `d2' / `tolse'
	}
	frame gf: local gl = estimate[`=`r0' + `k'']
	local d3 = reldif(e(ll), `gl')
	if `d3' / `tol' > `worst' local worst = `d3' / `tol'
	frame gf: drop __rn
	if `worst' > 1 {
		di as error "`id': FAIL (worst difference to goldfish = " %9.2e `worst' " times its tolerance)"
		exit 9
	}
	di as result "PASS `id': estimates, SEs, ll agree with goldfish (worst " %9.2e `worst' " of the tolerance), ll " %12.6f e(ll)
end
