* nwdynam postestimation against goldfish 1.6.12 (+ R package changepoint 2.3): reads the files of dev/dynam_post_crosscheck.R (global POSTDIR).
* Per case: estat ic (log likelihood, N, AIC, BIC), estat eventlik (log likelihood of every event, goldfish's intervalLogL; tolerance POSTTOL,
* default 1e-10 absolute), estat eventprob (goldfish's eventProbabilities, same tolerance), estat outliers and estat changepoints (exact indices).
* Globals: POSTONLY = only these cases, POSTNOPROB = 1 skips the (large) probability comparison.

capture program drop pimpd
program pimpd
	syntax using/ [, stringcols(string)]
	import delimited using `"`using'"', clear varnames(1) stringcols(`stringcols') asdouble
	quietly ds, has(type float)
	foreach v in `r(varlist)' {
		quietly recast double `v'
	}
end

* fit one case (same data preparation as the composition and exogenous cross-checks)
capture program drop prun
program prun
	args cid csub copts isdir exo
	tempfile ndfile
	pimpd using `"$POSTDIR/nodes_`cid'.csv"', stringcols(1)
	save `ndfile'
	nwclear
	pimpd using `"$POSTDIR/events_`cid'.csv"', stringcols(2 3)
	if "`isdir'" == "2" nwset sender receiver, twomode eventtime(time) name(net)
	else if "`isdir'" == "0" nwset sender receiver, eventtime(time) undirected name(net)
	else nwset sender receiver, eventtime(time) directed name(net)
	if "`isdir'" != "2" {
		_nwsyntax net, max(1)
		preserve
		use `ndfile', clear
		local missing ""
		forvalues i = 1/`=_N' {
			local lb = label[`i']
			mata: st_local("has", strofreal(any(`netobj'->get_nodenames() :== "`lb'")))
			if `has' == 0 local missing "`missing' `lb'"
		}
		restore
		foreach lb of local missing {
			nwaddnodes net, nodenames(`lb')
		}
	}
	capture confirm file `"$POSTDIR/changes_x_`cid'.csv"'
	if !_rc {
		preserve
		pimpd using `"$POSTDIR/changes_x_`cid'.csv"', stringcols(2)
		generate node2 = node
		nwset node node2 value, eventtime(time) selfloop name(chg_x)
		restore
	}
	capture confirm file `"$POSTDIR/comp_`cid'.csv"'
	if !_rc {
		preserve
		pimpd using `"$POSTDIR/comp_`cid'.csv"', stringcols(2)
		generate node2 = node
		nwset node node2 value, eventtime(time) selfloop name(comp)
		restore
	}
	_nwsyntax net, max(1)
	mata: st_local("nact", strofreal(length(`netobj'->get_nodenames())))
	foreach e of local exo {
		local nm = substr("`e'", 1, strpos("`e'", "=") - 1)
		local ty = substr("`e'", strpos("`e'", "=") + 1, .)
		local dr "directed"
		if substr("`ty'", -3, 3) == "und" local dr "undirected"
		if substr("`ty'", 1, 3) == "dyn" {
			preserve
			pimpd using `"$POSTDIR/exo_`cid'_`nm'.csv"', stringcols(2 3)
			nwset sender receiver increment, eventtime(time) `dr' name(`nm')
			restore
		}
		else {
			preserve
			pimpd using `"$POSTDIR/exo_`cid'_`nm'.csv"'
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
	merge 1:1 label using `ndfile', nogen
	rename label _nwnode
	capture confirm file `"$POSTDIR/avail_`cid'.csv"'
	if !_rc {
		capture frame drop oppdata
		frame create oppdata
		frame oppdata: pimpd using `"$POSTDIR/avail_`cid'.csv"'
	}
	nwdynam net, submodel(`csub') `copts'
end

* the integer list of a row vector as a string
capture program drop vecstr
program vecstr, rclass
	args M
	local s ""
	forvalues j = 1/`=colsof(`M')' {
		local s "`s' `=el(`M', 1, `j')'"
	}
	return local s = trim("`s'")
end

capture program drop pcheck
program pcheck, eclass
	args cid tol
	* the log likelihood at goldfish's own estimates: the per-event values are compared at the same parameter vector
	frame pres: generate long __rn = _n
	frame pres: quietly summarize __rn if case == "`cid'"
	frame pres: local r0 = r(min)
	frame pres: local r1 = r(max)
	local k = colsof(e(b))
	if `r1' - `r0' != `k' {
		di as error "`cid': goldfish has " `r1' - `r0' " coefficients, nwdynam `k'"
		exit 9
	}
	tempname B0 B1
	matrix `B0' = e(b)
	matrix `B1' = e(b)
	local dest 0
	forvalues j = 1/`k' {
		frame pres: local gb = estimate[`=`r0' + `j' - 1']
		matrix `B1'[1, `j'] = `gb'
		local dest = max(`dest', abs(`B0'[1, `j'] - `gb') / (1 + abs(`gb')))
	}
	frame pres: local gll = estimate[`r1']
	frame pres: drop __rn
	if `dest' > 1e-5 {
		di as error "`cid': estimates differ from goldfish (`dest')"
		exit 9
	}
	local ownll = e(ll)
	ereturn repost b = `B1', rename
	ereturn scalar ll = `gll'
	* -- estat ic
	frame post_ic {
		quietly summarize ll if case == "`cid'"
		local gll2 = r(mean)
		quietly summarize nevents if case == "`cid'"
		local gn = r(mean)
		quietly summarize aic if case == "`cid'"
		local gaic = r(mean)
		quietly summarize bic if case == "`cid'"
		local gbic = r(mean)
		quietly summarize avg if case == "`cid'"
		local gavg = r(mean)
		quietly summarize ndep if case == "`cid'"
		local gdep = r(mean)
	}
	quietly estat ic
	local d = max(reldif(r(ll), `gll2'), reldif(r(aic), `gaic'), reldif(r(bic), `gbic'), reldif(r(ll_per_event), `gavg'))
	if r(N) != `gn' local d = 1
	if `d' > 1e-10 {
		di as error "`cid': estat ic differs from goldfish (N " r(N) " vs `gn'; ll " r(ll) " vs `gll'; aic " r(aic) " vs `gaic'; bic " r(bic) " vs `gbic')"
		exit 9
	}
	* -- estat eventlik
	capture frame drop pev
	quietly estat eventlik, frame(pev)
	if r(N) != `gdep' {
		di as error "`cid': estat eventlik reports " r(N) " events, goldfish `gdep'"
		exit 9
	}
	local ncens = r(n_cens)
	local lltrail = r(ll_trailing)
	local llev = r(ll_events)
	frame pev: generate long __k = _n
	frame pllref: quietly keep if case == "`cid'"
	capture frame drop pref
	frame pllref: frame put if case == "`cid'" & type == 1, into(pref)
	frame pref: generate long __k = _n
	frame pref: rename ll llref
	frame pev: frlink 1:1 __k, frame(pref)
	frame pev: frget llref, from(pref)
	frame pev: generate double __d = abs(ll - llref)
	frame pev: quietly summarize __d
	local dll = r(max)
	frame pllref: quietly count if case == "`cid'" & type == 2
	local gcens = r(N)
	frame pllref: quietly summarize ll if case == "`cid'" & type == 2
	local gsumcens = cond(r(N) > 0, r(sum), 0)
	frame pev: capture confirm variable llcens
	local hasc = !_rc
	local sumcens = `lltrail'
	if `hasc' {
		frame pev: quietly summarize llcens
		local sumcens = `sumcens' + cond(r(N) > 0, r(sum), 0)
	}
	if `ncens' != `gcens' | abs(`sumcens' - `gsumcens') > `tol' * (1 + abs(`gsumcens')) {
		di as error "`cid': censored intervals: ours `ncens' sum `sumcens', goldfish `gcens' sum `gsumcens'"
		exit 9
	}
	if `dll' > `tol' {
		di as error "`cid': per-event log likelihood differs from goldfish by " %9.2e `dll'
		exit 9
	}
	capture frame drop pref
	* -- estat eventprob
	local dpr 0
	if "$POSTNOPROB" == "" {
		capture frame drop ppr
		quietly estat eventprob, frame(ppr)
		capture frame drop pprref
		frame create pprref
		frame pprref: pimpd using `"$POSTDIR/ref_pr_`cid'.csv"'
		capture frame pprref: confirm variable actor1
		if !_rc {
			frame pprref: generate str40 k1 = cond(actor1 < actor2, actor1, actor2)
			frame pprref: generate str40 k2 = cond(actor1 < actor2, actor2, actor1)
			frame ppr: generate str40 k1 = cond(actor1 < actor2, actor1, actor2)
			frame ppr: generate str40 k2 = cond(actor1 < actor2, actor2, actor1)
			frame ppr: rename actor1 __a1
			frame ppr: rename actor2 __a2
			frame pprref: rename actor1 __b1
			frame pprref: rename actor2 __b2
			frame pprref: rename value ref
			frame ppr: keep event k1 k2 p
			frame pprref: keep event k1 k2 ref
			frame ppr: frlink 1:1 event k1 k2, frame(pprref)
		}
		else {
			frame pprref: rename value ref
			frame ppr: capture rename lnrate p
			frame ppr: keep event actor p
			frame pprref: keep event actor ref
			frame ppr: frlink 1:1 event actor, frame(pprref)
		}
		frame pprref: local nref = _N
		frame ppr: frget ref, from(pprref)
		frame ppr: quietly count if !missing(ref)
		if r(N) != `nref' {
			di as error "`cid': estat eventprob has " r(N) " rows matching goldfish's " `nref'
			exit 9
		}
		frame ppr: quietly count if !missing(p) & missing(ref)
		local nextra = r(N)
		frame ppr: generate double __d = abs(p - ref)
		frame ppr: quietly summarize __d
		local dpr = r(max)
		if `dpr' > `tol' {
			di as error "`cid': event probabilities differ from goldfish by " %9.2e `dpr'
			exit 9
		}
		capture frame drop ppr
		capture frame drop pprref
	}
	* -- outliers: exact indices
	local nol 0
	frame pout: frame put if case == "`cid'", into(pcase)
	frame pcase {
		local nr = _N
		forvalues i = 1/`nr' {
			local om`i' = lower(method[`i'])
			local op`i' = parameter[`i']
			local ow`i' = window[`i']
			local oi`i' = trim(idx[`i'])
		}
	}
	capture frame drop pcase
	forvalues i = 1/`nr' {
		local wopt ""
		if `ow`i'' >= 0 local wopt "window(`ow`i'')"
		capture quietly estat outliers, method(`om`i'') parameter(`op`i'') `wopt'
		local rc = _rc
		if "`oi`i''" == "ERR" {
			if `rc' == 0 {
				di as error "`cid': estat outliers `om`i'' `op`i'' `wopt' ran, goldfish gave an error"
				exit 9
			}
			continue
		}
		if `rc' {
			di as error "`cid': estat outliers `om`i'' `op`i'' `wopt' failed r(`rc'), goldfish `oi`i''"
			exit 9
		}
		local ours ""
		if r(n_outliers) > 0 {
			tempname M
			matrix `M' = r(outliers)
			vecstr `M'
			local ours = r(s)
		}
		* method top: events with the same log likelihood (to the last bits) may be ranked either way; goldfish and we then pick different ones
		if "`ours'" != "`oi`i''" & "`om`i''" == "top" {
			local same 1
			local ln : word count `ours'
			local ng : word count `oi`i''
			if `ln' == `ng' {
				local lo ""
				local lg ""
				foreach j of local ours {
					frame pev: local vv = ll[`j']
					local lo "`lo' `vv'"
				}
				foreach j of local oi`i' {
					frame pev: local vv = ll[`j']
					local lg "`lg' `vv'"
				}
				forvalues a = 1/`ln' {
					local va : word `a' of `lo'
					local found 0
					forvalues b = 1/`ln' {
						local vb : word `b' of `lg'
						if abs(`va' - `vb') < 1e-12 & `found' == 0 local found = `b'
					}
					if `found' == 0 local same 0
				}
			}
			else local same 0
			if `same' local ours "`oi`i''"
		}
		if "`ours'" != "`oi`i''" {
			di as error "`cid': outliers `om`i'' `op`i'' `wopt': ours [`ours'] goldfish [`oi`i'']"
			exit 9
		}
		local nol = `nol' + 1
	}
	* -- change points: exact
	frame pcpt: frame put if case == "`cid'", into(pcase)
	frame pcase {
		local nr = _N
		forvalues i = 1/`nr' {
			local cm`i' = lower(moment[`i'])
			local cg`i' = lower(method[`i'])
			local cw`i' = window[`i']
			local ci`i' = trim(idx[`i'])
		}
	}
	capture frame drop pcase
	local ncp 0
	forvalues i = 1/`nr' {
		local wopt ""
		if `cw`i'' > 0 local wopt "window(`cw`i'')"
		capture quietly estat changepoints, moment(`cm`i'') method(`cg`i'') `wopt'
		local rc = _rc
		if "`ci`i''" == "ERR" {
			if `rc' == 0 {
				di as error "`cid': estat changepoints `cm`i'' `cg`i'' `wopt' ran, goldfish gave an error"
				exit 9
			}
			continue
		}
		if `rc' {
			di as error "`cid': estat changepoints `cm`i'' `cg`i'' `wopt' failed r(`rc'), goldfish `ci`i''"
			exit 9
		}
		tempname C
		matrix `C' = r(cpts)
		vecstr `C'
		if "`r(s)'" != "`ci`i''" {
			di as error "`cid': changepoints `cm`i'' `cg`i'' `wopt': ours [`r(s)'] goldfish [`ci`i'']"
			exit 9
		}
		local ncp = `ncp' + 1
	}
	capture frame drop pev
	di as result "PASS `cid': ic, eventlik (max diff " %8.1e `dll' "), eventprob (max diff " %8.1e `dpr' "), `nol' outlier and `ncp' change point settings exact; censored intervals `ncens'"
end

capture program drop prunall
program prunall
	capture frame drop post_ic
	capture frame drop pres
	capture frame drop pllref
	capture frame drop pout
	capture frame drop pcpt
	frame create post_ic
	frame post_ic: pimpd using `"$POSTDIR/post_ic.csv"', stringcols(1)
	frame create pout
	frame pout: pimpd using `"$POSTDIR/post_outliers.csv"', stringcols(1 2 5)
	frame create pcpt
	frame pcpt: pimpd using `"$POSTDIR/post_cpts.csv"', stringcols(1 2 3 5)
	* the reference log likelihoods of all cases in one frame
	frame create pllref
	capture frame drop pres
	frame create pres
	frame pres: pimpd using `"$POSTDIR/results.csv"', stringcols(1 2)
	preserve
	pimpd using `"$POSTDIR/cases.csv"', stringcols(1 2 3 7)
	local nc = _N
	forvalues i = 1/`nc' {
		local id`i' = id[`i']
		local sub`i' = sub[`i']
		local opts`i' = opts[`i']
		local dir`i' = directed[`i']
		local exo`i' = exo[`i']
	}
	restore
	local tol : env POSTTOL
	if "`tol'" == "" local tol 1e-10
	local nfail 0
	forvalues i = 1/`nc' {
		if "$POSTONLY" != "" & !strpos(" $POSTONLY ", " `id`i'' ") continue
		capture noisily prun `id`i'' `sub`i'' "`opts`i''" `dir`i'' "`exo`i''"
		if _rc {
			di as error "`id`i'': nwdynam failed, r(" _rc ")"
			local nfail = `nfail' + 1
			continue
		}
		frame pllref: clear
		frame pllref {
			quietly import delimited using `"$POSTDIR/ref_ll_`id`i''.csv"', clear varnames(1) asdouble
			quietly generate str20 case = "`id`i''"
		}
		capture noisily pcheck `id`i'' `tol'
		if _rc local nfail = `nfail' + 1
	}
	if `nfail' > 0 {
		di as error "`nfail' case(s) differ from goldfish"
		exit 9
	}
	di as result "ALL OK: `nc' cases agree with goldfish"
end
