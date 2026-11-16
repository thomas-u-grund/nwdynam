* nwdynam effects on other networks, mixed closures, tertiusDiff, tie() of changing networks and ignoreRep against goldfish: reads the files of
* dev/dynam_ex_crosscheck.R (global WINDIR) and compares estimates, standard errors and log likelihood (tolerance WINTOL, default 1e-6; choice_coordination:
* estimates 5e-4, standard errors 10%, as their standard errors come from a numerical derivative).
* The caller sets the paths (adopath of the core and of this package, a scratch PLUS directory) and the global WINDIR.
* Globals: WINSEC = the sections to run (a b c d e f g; default all with a cases file), WINONLY = only these cases.

* import delimited stores 2-digit decimals as float; goldfish reads doubles (window expiry times are sensitive to that)
capture program drop impd
program impd
	syntax using/ [, stringcols(string)]
	import delimited using `"`using'"', clear varnames(1) stringcols(`stringcols')
	quietly ds, has(type float)
	foreach v in `r(varlist)' {
		quietly recast double `v'
	}
end

capture program drop runcase
program runcase
	args id sub opts kind exo
	tempfile nodes
	impd using `"$WINDIR/nodes_`id'.csv"', stringcols(1)
	save `nodes'
	local nact = _N
	nwclear
	impd using `"$WINDIR/events_`id'.csv"', stringcols(2 3)
	local vv ""
	capture confirm variable increment
	if !_rc local vv "increment"
	if "`kind'" == "und" nwset sender receiver `vv', eventtime(time) undirected name(net)
	else nwset sender receiver `vv', eventtime(time) directed name(net)
	foreach e of local exo {
		local nm = substr("`e'", 1, strpos("`e'", "=") - 1)
		local ty = substr("`e'", strpos("`e'", "=") + 1, .)
		local dr "directed"
		if substr("`ty'", -3, 3) == "und" local dr "undirected"
		if substr("`ty'", 1, 3) == "dyn" {
			preserve
			impd using `"$WINDIR/exo_`id'_`nm'.csv"', stringcols(2 3)
			nwset sender receiver increment, eventtime(time) `dr' name(`nm')
			restore
		}
		else {
			preserve
			impd using `"$WINDIR/exo_`id'_`nm'.csv"'
			mata: __exm = J(`nact', `nact', 0)
			mata: __exd = st_data(., "sender receiver value")
			mata: for (__exi = 1; __exi <= rows(__exd); __exi++) __exm[__exd[__exi, 1], __exd[__exi, 2]] = __exd[__exi, 3]
			mata: st_matrix("__exmat", __exm)
			mata: mata drop __exm __exd __exi
			restore
			preserve
			nwset, mat(__exmat) name(`nm') `dr' selfloop
			restore
		}
	}
	rename _nwnode label
	merge 1:1 label using `nodes', keep(match master) nogen
	rename label _nwnode
	nwdynam net, submodel(`sub') `opts'
end

capture program drop chkcase
program chkcase
	args id tolb tolse
	local tol : env WINTOL
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
		di as error "`id': FAIL (worst difference to goldfish = " %9.2e `worst' " times its tolerance; ll " %12.6f e(ll) " vs " %12.6f `gl' ")"
		exit 9
	}
	di as result "PASS `id': estimates, SEs, ll agree with goldfish (worst " %9.2e `worst' " of the tolerance), ll " %12.6f e(ll)
end

capture program drop runsec
program runsec
	args sec
	preserve
	import delimited using `"$WINDIR/results_`sec'.csv"', clear varnames(1)
	tempfile res
	save `res'
	capture frame drop gf
	frame create gf
	frame gf: use `res', clear
	import delimited using `"$WINDIR/cases_`sec'.csv"', clear varnames(1) stringcols(1 2 3 4 5)
	local nc = _N
	forvalues i = 1/`nc' {
		local id`i' = id[`i']
		local sub`i' = sub[`i']
		local opts`i' = opts[`i']
		local kind`i' = kind[`i']
		local exo`i' = exo[`i']
		local tb`i' = tolb[`i']
		local ts`i' = tolse[`i']
	}
	restore
	local nfail 0
	forvalues i = 1/`nc' {
		if "$WINONLY" != "" & !strpos(" $WINONLY ", " `id`i'' ") continue
		capture noisily runcase `id`i'' `sub`i'' `"`opts`i''"' `kind`i'' "`exo`i''"
		if _rc {
			di as error "`id`i'': nwdynam failed, r(" _rc ")"
			local nfail = `nfail' + 1
			continue
		}
		capture noisily chkcase `id`i'' `tb`i'' `ts`i''
		if _rc local nfail = `nfail' + 1
	}
	if `nfail' > 0 {
		di as error "section `sec': `nfail' case(s) differ from goldfish"
		exit 9
	}
	di as result "section `sec': ALL OK, `nc' cases agree with goldfish"
end

capture program drop runall
program runall
	local secs "$WINSEC"
	if "`secs'" == "" local secs "a b c d e f g"
	local nfail 0
	foreach s of local secs {
		capture confirm file `"$WINDIR/cases_`s'.csv"'
		if _rc continue
		capture noisily runsec `s'
		if _rc local nfail = `nfail' + 1
	}
	if `nfail' > 0 {
		di as error "`nfail' section(s) with differences"
		exit 9
	}
	di as result "ALL OK"
end
