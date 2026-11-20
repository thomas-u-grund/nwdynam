*! version 1.9.1  18nov2026 nwdynam_estat
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
* The change point methods of estat changepoints are a port of the R package 'changepoint' (Rebecca Killick, Kaylea Haynes, Idris Eckley, Paul Fearnhead, Jamie Lee; GPL), see NOTICE.md.
* nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

* Postestimation of nwdynam: estat ic, estat eventlik, estat eventprob, estat outliers, estat changepoints, and (extensions, not in goldfish) estat choiceplot and estat rateplot.
* What goldfish offers after estimation and nothing beyond it: logLik / AIC / BIC, the estimation options returnIntervalLogL
* (estat eventlik) and returnEventProbabilities (estat eventprob), examineOutliers (estat outliers), examineChangepoints
* (estat changepoints).

capture program drop nwdynam_estat
program define nwdynam_estat, rclass
	gettoken subcmd rest : 0, parse(" ,")
	local subcmd = lower(trim("`subcmd'"))
	if `"`e(cmd)'"' != "nwdynam" {
		di as err "last nwdynam estimates not found"
		exit 301
	}
	_nwdynam_libcheck
	if "`subcmd'" == "ic" {
		nwdynam_estat_ic `rest'
		return add
	}
	else if "`subcmd'" == "eventlik" {
		nwdynam_estat_eventlik `rest'
		return add
	}
	else if "`subcmd'" == "eventprob" {
		nwdynam_estat_eventprob `rest'
		return add
	}
	else if "`subcmd'" == "outliers" {
		nwdynam_estat_outliers `rest'
		return add
	}
	else if "`subcmd'" == "changepoints" {
		nwdynam_estat_changepoints `rest'
		return add
	}
	else if "`subcmd'" == "choiceplot" {
		nwdynam_estat_choiceplot `rest'
		return add
	}
	else if "`subcmd'" == "rateplot" {
		nwdynam_estat_rateplot `rest'
		return add
	}
	else {
		// estat vce, estat summarize ... : Stata's defaults; anything else is refused with r(321)
		estat_default `0'
	}
end

* the fit that estat works on must still be in memory and be the one of the stored results
capture program drop _nwdynam_estat_check
program define _nwdynam_estat_check
	capture confirm frame __nwdynam_fit
	local ok = (_rc == 0)
	if `ok' {
		mata: st_local("__nwdynam_fid", dynam_fit_id())
		if `"`e(fitid)'"' == "" | `"`__nwdynam_fid'"' != `"`e(fitid)'"' local ok = 0
	}
	if !`ok' {
		di as err "{p 0 0 2}no fitted nwdynam model available for this command (the data of the fit is kept in memory after {bf:nwdynam} and is gone after {bf:clear all}, a later fit, or when the estimates come from {bf:estimates restore} of another fit). Run {bf:nwdynam} again.{p_end}"
		exit 498
	}
end

* a frame name for the result: must not exist unless replace is given
capture program drop _nwdynam_estat_frame
program define _nwdynam_estat_frame
	args fname replace
	capture confirm name `fname'
	if _rc {
		di as err "frame(): {bf:`fname'} is not a valid frame name"
		exit 198
	}
	capture confirm frame `fname'
	if !_rc {
		if "`replace'" == "" {
			di as err "{p 0 0 2}frame {bf:`fname'} already exists; use option {bf:replace}.{p_end}"
			exit 110
		}
		if "`fname'" == "`c(frame)'" {
			di as err "frame(): {bf:`fname'} is the current frame"
			exit 198
		}
		frame drop `fname'
	}
end

* ---------------------------------------------------------------------------------------------------------------------
capture program drop nwdynam_estat_ic
program define nwdynam_estat_ic, rclass
	syntax
	_nwdynam_estat_check
	local fr = c(frame)
	capture noisily mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(0, 0), "%18.0f"))
	if _rc {
		local rc = _rc
		capture frame change `fr'
		exit `rc'
	}
	capture frame change `fr'
	// goldfish: nobs = nEvents = every record of the preprocessing, in the rate model with intercept also its right-censored intervals
	local N = `__nwdynam_nrec'
	local df = colsof(e(b))
	local ll = e(ll)
	local aic = -2 * `ll' + 2 * `df'
	local bic = -2 * `ll' + ln(`N') * `df'
	local avg = `ll' / `N'
	tempname M
	matrix `M' = (`N', `ll', `avg', `df', `aic', `bic')
	matrix colnames `M' = N ll ll/N df AIC BIC
	matrix rownames `M' = `e(depvar)'
	di
	di as txt "Akaike's information criterion and Bayesian information criterion"
	di
	di as txt "{hline 13}{c TT}{hline 66}"
	di as txt %12s "Model" " {c |}" %9s "N" %13s "ll" %11s "ll/N" %6s "df" %11s "AIC" %11s "BIC"
	di as txt "{hline 13}{c +}{hline 66}"
	di as txt %12s abbrev("`e(depvar)'", 12) " {c |}" as res %9.0g `N' %13.7g `ll' %11.6g `avg' %6.0f `df' %11.7g `aic' %11.7g `bic'
	di as txt "{hline 13}{c BT}{hline 66}"
	if `N' != e(N) {
		di as txt "{p 0 6 2}Note: N counts every interval of the log likelihood, as goldfish does: " as res e(N) as txt " events and " as res `N' - e(N) as txt " right-censored intervals (periods ended by a covariate change or a window expiry). BIC uses this N.{p_end}"
	}
	else di as txt "Note: BIC uses N = number of events."
	return scalar ll = `ll'
	return scalar ll_per_event = `avg'
	return scalar N = `N'
	return scalar N_events = e(N)
	return scalar df = `df'
	return scalar aic = `aic'
	return scalar bic = `bic'
	return matrix S = `M'
end

* ---------------------------------------------------------------------------------------------------------------------
capture program drop nwdynam_estat_eventlik
program define nwdynam_estat_eventlik, rclass
	syntax [, FRAme(name) REPLACE]
	_nwdynam_estat_check
	if "`frame'" != "" _nwdynam_estat_frame `frame' `replace'
	local fr = c(frame)
	tempname tab
	frame create `tab'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(1, 0), "%18.0f"))
		mata: dynam_post_table("`tab'", `__nwdynam_nrec')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `tab'
		exit `rc'
	}
	capture frame change `fr'
	frame `tab' {
		quietly summarize ll
		local N = r(N)
		local sum = r(sum)
		local mean = r(mean)
		local min = r(min)
		local max = r(max)
	}
	di
	di as txt "Log likelihood of every event (goldfish: returnIntervalLogL), submodel " as res "`e(submodel)'" as txt cond("`e(effects)'" != "", "", "")
	di as txt "  Events: " as res `N' as txt "   sum: " as res %12.6f `sum' as txt "   mean: " as res %9.5f `mean' as txt "   min: " as res %9.5f `min' as txt "   max: " as res %9.5f `max'
	if `__nwdynam_hascens' {
		di as txt "{p 2 2 2}Rate sub-model with intercept: " as res `__nwdynam_ncens' as txt " right-censored intervals (periods ended by a covariate change or a window expiry) are part of the log likelihood as in goldfish; variable {bf:llcens} holds those that end before an event, " as res %9.5f `__nwdynam_lltrail' as txt " lie after the last event.{p_end}"
	}
	return scalar N = `N'
	return scalar ll_events = `sum'
	return scalar ll = e(ll)
	return scalar ll_mean = `mean'
	return scalar ll_min = `min'
	return scalar ll_max = `max'
	return scalar n_cens = `__nwdynam_ncens'
	return scalar ll_trailing = `__nwdynam_lltrail'
	if "`frame'" != "" {
		frame copy `tab' `frame'
		di as txt "Table of the events in frame " as res "`frame'" as txt ": event, origevent, time, sender, receiver, ll" cond(`__nwdynam_hascens', ", llcens", "") "."
		return local frame "`frame'"
	}
	capture frame drop `tab'
end

* ---------------------------------------------------------------------------------------------------------------------
capture program drop nwdynam_estat_eventprob
program define nwdynam_estat_eventprob, rclass
	syntax , FRAme(name) [REPLACE MAXRows(real 200000000)]
	_nwdynam_estat_check
	_nwdynam_estat_frame `frame' `replace'
	local fr = c(frame)
	tempname tab
	frame create `tab'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(2, `maxrows'), "%18.0f"))
		mata: dynam_post_prob("`tab'", `__nwdynam_nrec')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `tab'
		exit `rc'
	}
	capture frame change `fr'
	frame `tab': local nrows = _N
	frame copy `tab' `frame'
	capture frame drop `tab'
	di
	di as txt "Model probabilities of every event (goldfish: returnEventProbabilities), submodel " as res "`e(submodel)'" as txt ": " as res `nrows' as txt " rows in frame " as res "`frame'"
	if "`e(submodel)'" == "choice" di as txt "{p 2 2 2}One row per event and candidate receiver: {bf:p} = probability that the actor is the receiver; {bf:observed} = the realised receiver.{p_end}"
	else if "`e(submodel)'" == "choice_coordination" di as txt "{p 2 2 2}One row per event and pair of actors (actor1 < actor2): {bf:p} = goldfish's matrix entry, the probability of the pair times two (the sum over the pairs of an event is 1); {bf:observed} = the realised pair.{p_end}"
	else if "`e(effects)'" != "" & strpos("`e(effects)'", "Intercept") > 0 di as txt "{p 2 2 2}One row per event and actor: {bf:lnrate} = the actor's linear predictor minus the exposure of the period that ends at the event (goldfish returns this quantity for the rate model with intercept; it is not a probability); {bf:observed} = the sender.{p_end}"
	else di as txt "{p 2 2 2}One row per event and actor: {bf:p} = probability that the actor is the sender; {bf:observed} = the realised sender.{p_end}"
	return scalar N_rows = `nrows'
	return local frame "`frame'"
end

* ---------------------------------------------------------------------------------------------------------------------
capture program drop nwdynam_estat_outliers
program define nwdynam_estat_outliers, rclass
	syntax [, METhod(string) PARameter(real 3) WINdow(string) FRAme(name) REPLACE PLOT NAme(string)]
	_nwdynam_estat_check
	if "`method'" == "" local method hampel
	local method = lower("`method'")
	if !inlist("`method'", "hampel", "iqr", "top") {
		di as err "method() must be hampel, iqr or top"
		exit 198
	}
	local hasw = 0
	if `"`window'"' != "" {
		capture confirm number `window'
		if _rc | real("`window'") < 0 {
			di as err "window() must be a number of events, 0 or more"
			exit 198
		}
		local hasw = 1
	}
	if "`frame'" != "" _nwdynam_estat_frame `frame' `replace'
	local fr = c(frame)
	tempname tab
	frame create `tab'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(1, 0), "%18.0f"))
		mata: dynam_post_table("`tab'", `__nwdynam_nrec')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `tab'
		exit `rc'
	}
	capture frame change `fr'
	frame `tab': local N = _N
	local kind = cond("`method'" == "hampel", 1, cond("`method'" == "iqr", 2, 3))
	local wnum = cond(`hasw', real("`window'"), -1)
	if "`method'" == "top" {
		if `parameter' > `N' {
			capture frame drop `tab'
			di as err "{p 0 0 2}parameter() for method(top) is the number of events to flag and cannot exceed the number of events (" `N' ").{p_end}"
			exit 198
		}
	}
	mata: st_numscalar("__nwdynam_serr", dynam_post_series("`tab'", `kind', 0, 0, 0, `parameter', `wnum'))
	local err = __nwdynam_serr
	scalar drop __nwdynam_serr
	if `err' {
		capture frame drop `tab'
		di as err "{p 0 0 2}window() is too large for the number of events (" `N' "): 2*window + 1 must not exceed the number of events (goldfish stops with an error here as well).{p_end}"
		exit 198
	}
	local nout = `__nwdynam_nidx'
	di
	di as txt "Outliers of the log likelihood of the events, method " as res "`method'" as txt ", parameter " as res `parameter' as txt cond("`method'" == "hampel", ", window " + cond(`hasw', "`window'", string(`N'/2 - 1) + " (default, as goldfish)"), "")
	if "`method'" == "hampel" & !`hasw' {
		di as txt "{p 2 2 2}Note: with the default window the rule tests only the middle event(s) of the log; give window() for a useful rolling rule.{p_end}"
	}
	if `nout' == 0 {
		di as txt "No outliers found."
	}
	else {
		frame `tab' {
			generate strL label = sender + "-" + receiver if flag
			list event time sender receiver ll if flag, noobs sep(0) abbreviate(12)
		}
	}
	if `nout' > 0 {
		tempname idx
		matrix `idx' = __nwdynam_idx
		return matrix outliers = `idx'
	}
	capture matrix drop __nwdynam_idx
	return scalar n_outliers = `nout'
	return scalar N = `N'
	return scalar parameter = `parameter'
	if "`method'" == "hampel" return scalar window = cond(`hasw', real("`window'"), `N'/2 - 1)
	return local method "`method'"
	if "`plot'" != "" & `nout' >= 0 {
		frame `tab' {
			capture generate strL label = sender + "-" + receiver if flag
			local gname = cond("`name'" != "", "name(`name', replace)", "")
			twoway (line ll time, lcolor(black)) (scatter ll time if !flag, mcolor(black) msize(small)) ///
				(scatter ll time if flag, mcolor(red) msize(small) mlabel(label) mlabcolor(red) mlabsize(vsmall) mlabposition(6)), ///
				legend(off) xtitle("") ytitle("Interval log likelihood") `gname'
		}
	}
	if "`frame'" != "" {
		frame `tab': capture drop label
		frame `tab': rename flag outlier
		frame copy `tab' `frame'
		di as txt "Table of the events in frame " as res "`frame'" as txt " (variable outlier = 1 for the flagged events)."
		return local frame "`frame'"
	}
	capture frame drop `tab'
end

* ---------------------------------------------------------------------------------------------------------------------
capture program drop nwdynam_estat_changepoints
program define nwdynam_estat_changepoints, rclass
	syntax [, MOMent(string) METhod(string) WINdow(string) FRAme(name) REPLACE PLOT NAme(string)]
	_nwdynam_estat_check
	if "`moment'" == "" local moment mean
	local moment = lower("`moment'")
	if !inlist("`moment'", "mean", "variance") {
		di as err "moment() must be mean or variance"
		exit 198
	}
	if "`method'" == "" local method pelt
	local method = lower("`method'")
	if !inlist("`method'", "pelt", "amoc", "binseg") {
		di as err "method() must be pelt, amoc or binseg"
		exit 198
	}
	local hasw = 0
	if `"`window'"' != "" {
		capture confirm integer number `window'
		if _rc | real("`window'") < 1 {
			di as err "window() must be a whole number of events, 1 or more (the minimum segment length)"
			exit 198
		}
		local hasw = 1
	}
	if "`frame'" != "" _nwdynam_estat_frame `frame' `replace'
	local fr = c(frame)
	tempname tab
	frame create `tab'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(1, 0), "%18.0f"))
		mata: dynam_post_table("`tab'", `__nwdynam_nrec')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `tab'
		exit `rc'
	}
	capture frame change `fr'
	frame `tab': local N = _N
	local mnum = cond("`method'" == "pelt", 1, cond("`method'" == "amoc", 2, 3))
	local momnum = cond("`moment'" == "mean", 1, 2)
	local wnum = cond(`hasw', real("`window'"), -1)
	local __nwdynam_window = `wnum'
	mata: st_numscalar("__nwdynam_serr", dynam_post_series("`tab'", 4, `mnum', `momnum', max((`wnum', 1)), 0, `wnum'))
	local err = __nwdynam_serr
	scalar drop __nwdynam_serr
	if `err' {
		capture frame drop `tab'
		if `err' == 2 di as err "{p 0 0 2}too few events for a change point analysis (" `N' "; at least 2, for moment(variance) at least 4).{p_end}"
		else if `err' == 1 di as err "{p 0 0 2}window() is the minimum segment length and is too large for the number of events (" `N' "): the log needs at least 2*window events.{p_end}"
		else di as err "{p 0 0 2}method(binseg) looks for up to 5 change points, which needs more events or a smaller window() than this log has (" `N' " events).{p_end}"
		exit 198
	}
	local win = `__nwdynam_window'
	if "`moment'" == "variance" & `win' < 2 local win = 2
	local nidx = `__nwdynam_nidx'
	tempname idx
	matrix `idx' = __nwdynam_idx
	capture matrix drop __nwdynam_idx
	// the last entry is the end of the data (the change point list of the changepoint package ends with n); it is not a change
	local ncp = `nidx' - 1
	di
	di as txt "Change points of the log likelihood of the events, " as res "`method'" as txt " method, change in the " as res "`moment'" as txt ", minimum segment length (window) " as res `win'
	if `ncp' <= 0 {
		di as txt "No regime changes found."
	}
	else {
		frame `tab' {
			generate byte cpt = 0
			forvalues j = 1/`ncp' {
				replace cpt = 1 in `=`idx'[1, `j']'
			}
			list event time sender receiver ll if cpt, noobs sep(0) abbreviate(12)
		}
		di as txt "{p 2 2 2}The change points are the events that end a regime: the log likelihood of the following events belongs to the next regime.{p_end}"
	}
	// the time of every change point for the graph, read before `idx' is moved into r() (return matrix removes it)
	local xl ""
	if `ncp' > 0 {
		forvalues j = 1/`ncp' {
			local __nwdynam_row = `idx'[1, `j']
			frame `tab': local __nwdynam_tx = time[`__nwdynam_row']
			local xl "`xl' `__nwdynam_tx'"
		}
	}
	if `nidx' > 0 {
		return matrix cpts = `idx'
	}
	return scalar n_cpts = max(`ncp', 0)
	return scalar N = `N'
	return scalar window = `win'
	return local method "`method'"
	return local moment "`moment'"
	if "`plot'" != "" {
		frame `tab' {
			local gname = cond("`name'" != "", "name(`name', replace)", "")
			twoway (line ll time, lcolor(black)) (scatter ll time, mcolor(black) msize(small)), ///
				legend(off) xtitle("") ytitle("Interval log likelihood") `=cond("`xl'" != "", "xline(`xl', lcolor(red))", "")' `gname'
		}
	}
	if "`frame'" != "" {
		frame `tab' {
			capture confirm variable cpt
			if _rc generate byte cpt = 0
		}
		frame copy `tab' `frame'
		di as txt "Table of the events in frame " as res "`frame'" as txt " (variable cpt = 1 for the events that end a regime)."
		return local frame "`frame'"
	}
	capture frame drop `tab'
end

* ---------------------------------------------------------------------------------------------------------------------
* Extension (goldfish has no such plot): where does the realised receiver (pair) rank among the risk set of its event by the fitted probability
capture program drop nwdynam_estat_choiceplot
program define nwdynam_estat_choiceplot, rclass
	syntax [, TOPk(numlist integer >=1 sort) MAXRANk(integer 20) FRAme(name) REPLACE NOPLOT NAme(string) MAXRows(real 200000000) *]
	local twopt `options'
	_nwdynam_estat_check
	if !inlist("`e(submodel)'", "choice", "choice_coordination") {
		di as err "{p 0 0 2}estat choiceplot needs the sub-model {bf:choice} or {bf:choice_coordination} (the last fit is {bf:`e(submodel)'}); use {bf:estat rateplot} after the rate model.{p_end}"
		exit 198
	}
	if `maxrank' < 1 {
		di as err "maxrank() must be 1 or more"
		exit 198
	}
	if "`topk'" == "" local topk "1 3 5 10"
	if "`frame'" != "" _nwdynam_estat_frame `frame' `replace'
	local K = `maxrank'
	foreach k of local topk {
		if `k' > `K' local K = `k'
	}
	local fr = c(frame)
	tempname rk curve
	frame create `rk'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(2, `maxrows'), "%18.0f"))
		mata: dynam_post_rank("`rk'", "`curve'", `__nwdynam_nrec', `K')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `rk'
		exit `rc'
	}
	capture frame change `fr'
	frame `rk' {
		quietly count if rank < .
		local N = r(N)
		quietly summarize rank, detail
		local mrank = r(mean)
		local medrank = r(p50)
		quietly generate double __u = (m + 1) / 2 if rank < .
		quietly summarize __u
		local urank = r(mean)
		quietly summarize m
		local mrisk = r(mean)
		local maxm = r(max)
		quietly generate double __ratio = rank / m
		quietly summarize __ratio
		local mrel = r(mean)
		drop __u __ratio
	}
	if `N' == 0 {
		capture frame drop `rk'
		di as err "no event with a defined risk set"
		exit 198
	}
	di
	di as txt "Rank of the realised " cond("`e(submodel)'" == "choice", "receiver", "pair") " among the risk set by the fitted probability (rank 1 = most probable), " as res `N' as txt " events"
	di as txt "{p 2 2 2}Mean rank " as res %9.3f `mrank' as txt "   median rank " as res %7.1f `medrank' as txt "   expected mean rank if all risk-set members were equally probable " as res %9.3f `urank' as txt "   mean risk-set size " as res %9.2f `mrisk' as txt ".{p_end}"
	di
	di as txt "{hline 13}{c TT}{hline 47}"
	di as txt %12s "Top k" " {c |}" %15s "observed" %15s "uniform" %15s "model"
	di as txt "{hline 13}{c +}{hline 47}"
	tempname T
	local nt : word count `topk'
	matrix `T' = J(`nt', 4, .)
	local row 0
	foreach k of local topk {
		local ++row
		matrix `T'[`row', 1] = `k'
		matrix `T'[`row', 2] = `curve'[`k', 2]
		matrix `T'[`row', 3] = `curve'[`k', 3]
		matrix `T'[`row', 4] = `curve'[`k', 4]
		di as txt %12.0f `k' " {c |}" as res %15.4f `curve'[`k', 2] %15.4f `curve'[`k', 3] %15.4f `curve'[`k', 4]
	}
	di as txt "{hline 13}{c BT}{hline 47}"
	matrix colnames `T' = k observed uniform model
	if "`noplot'" == "" {
		local kp = min(`maxrank', max(`maxm', 1))
		tempname cm
		frame create `cm'
		frame `cm' {
			quietly set obs `kp'
			quietly generate double k = _n
			quietly generate double observed = .
			quietly generate double uniform = .
			quietly generate double model = .
			forvalues j = 1/`kp' {
				quietly replace observed = `curve'[`j', 2] in `j'
				quietly replace uniform = `curve'[`j', 3] in `j'
				quietly replace model = `curve'[`j', 4] in `j'
			}
			local gname = cond("`name'" != "", "name(`name', replace)", "")
			local what = cond("`e(submodel)'" == "choice", "receiver", "pair")
			twoway (line uniform k, lcolor(gs10) lpattern(dash) lwidth(medthick)) ///
				(line model k, lcolor(cranberry) lpattern(shortdash) lwidth(medthick)) ///
				(connected observed k, lcolor(navy) mcolor(navy) msize(small) lwidth(medthick)), ///
				ylabel(0(.2)1, angle(0) grid glcolor(gs14)) xlabel(, grid glcolor(gs14) format(%9.0f)) ///
				xtitle("Rank k of the realised `what'") ytitle("Share of events with rank {&le} k") ///
				legend(order(3 "Observed" 2 "Model expectation" 1 "Uniform (all equally probable)") rows(1) size(vsmall) symxsize(5)) ///
				title("Fitted probability of the realised `what': rank in the risk set", size(medsmall)) ///
				graphregion(color(white)) plotregion(margin(small)) `gname' `twopt'
		}
		capture frame drop `cm'
	}
	return scalar N = `N'
	return scalar mean_rank = `mrank'
	return scalar median_rank = `medrank'
	return scalar mean_rank_uniform = `urank'
	return scalar mean_rank_relative = `mrel'
	return scalar mean_riskset = `mrisk'
	return matrix topk = `T'
	tempname C
	matrix `C' = `curve'
	matrix colnames `C' = k observed uniform model
	return matrix curve = `C'
	if "`frame'" != "" {
		frame copy `rk' `frame'
		di as txt "Table of the events in frame " as res "`frame'" as txt ": event, time, m (size of the risk set), rank, pobs (fitted probability of the realised " cond("`e(submodel)'" == "choice", "receiver", "pair") ")."
		return local frame "`frame'"
	}
	capture frame drop `rk'
end

* ---------------------------------------------------------------------------------------------------------------------
* Extension (goldfish has no such plot): observed and expected number of events per actor as sender, and per period
capture program drop nwdynam_estat_rateplot
program define nwdynam_estat_rateplot, rclass
	syntax [, OVer(string) PERiods(integer 10) TOPActors(integer 5) LABel(integer 5) FRAme(name) REPLACE NOPLOT NAme(string) MAXRows(real 200000000) *]
	local twopt `options'
	_nwdynam_estat_check
	if "`e(submodel)'" != "rate" {
		di as err "{p 0 0 2}estat rateplot needs the sub-model {bf:rate} (the last fit is {bf:`e(submodel)'}); use {bf:estat choiceplot} after the choice model.{p_end}"
		exit 198
	}
	if "`over'" == "" local over actor
	local over = lower("`over'")
	if !inlist("`over'", "actor", "time") {
		di as err "over() must be actor or time"
		exit 198
	}
	if `periods' < 1 {
		di as err "periods() must be 1 or more"
		exit 198
	}
	if `topactors' < 1 | `label' < 0 {
		di as err "topactors() must be 1 or more and label() 0 or more"
		exit 198
	}
	if "`frame'" != "" _nwdynam_estat_frame `frame' `replace'
	local fr = c(frame)
	tempname act per
	frame create `act'
	frame create `per'
	capture noisily {
		mata: st_local("__nwdynam_nrec", strofreal(dynam_post_run(2, `maxrows'), "%18.0f"))
		mata: dynam_post_rate("`act'", "`per'", `__nwdynam_nrec', `periods')
	}
	if _rc {
		local rc = _rc
		capture frame change `fr'
		capture frame drop `act'
		capture frame drop `per'
		exit `rc'
	}
	capture frame change `fr'
	local hasint = (strpos("`e(effects)'", "Intercept") > 0)
	frame `act' {
		quietly count
		local n = r(N)
		quietly summarize observed
		local tobs = r(sum)
		quietly correlate observed expected
		local corr = cond(r(N) > 1 & r(rho) < ., r(rho), .)
		quietly generate double __d = abs(observed - expected)
		quietly generate double __sd = (observed - expected) / sqrt(max(expected, 1e-12))
		quietly summarize __sd
		quietly generate double __chi = __sd^2
		quietly summarize __chi if expected > 0
		local chi = r(sum)
		local nchi = r(N)
	}
	di
	di as txt "Observed and expected number of events per actor as sender, " as res `n' as txt " actors, " as res trim(string(`tobs', "%12.0fc")) as txt " events"
	di as txt "{p 2 2 2}Correlation of observed and expected counts " as res %7.4f `corr' as txt "   Pearson statistic sum (obs-exp)^2/exp " as res %10.2f `chi' as txt " over " as res `nchi' as txt " actors (descriptive only, not a test).{p_end}"
	if `hasint' di as txt "{p 2 2 2}Rate model with intercept: the expected count of an actor is the sum over the events of its share of the total hazard, exp(lnrate) divided by the sum over the actors; the expected counts add up to the observed number of events, so the plot checks who acts, not how many events occur.{p_end}"
	if "`over'" == "actor" {
		frame `act' {
			gsort -expected
			local kk = min(`label', `n')
			format observed %9.0f
			format expected %9.2f
			list actor observed expected in 1/`=min(10, `n')', noobs sep(0) abbreviate(14)
		}
		if "`noplot'" == "" {
			frame `act' {
				quietly generate strL __lab = ""
				gsort -__d
				if `kk' > 0 quietly replace __lab = actor in 1/`kk'
				quietly summarize observed
				local mx = r(max)
				quietly summarize expected
				local mx = max(`mx', r(max))
				local gname = cond("`name'" != "", "name(`name', replace)", "")
				twoway (function y = x, range(0 `mx') lcolor(gs10) lpattern(dash) lwidth(medthick)) ///
					(scatter observed expected, mcolor(navy%60) msize(small) mlabel(__lab) mlabcolor(navy) mlabsize(vsmall) mlabposition(1)), ///
					aspectratio(1) xtitle("Expected number of events as sender") ytitle("Observed number of events as sender") ///
					xlabel(, grid glcolor(gs14) format(%9.0f)) ylabel(, angle(0) grid glcolor(gs14) format(%9.0f)) legend(off) ///
					title("Sender activity: observed against expected under the model", size(medsmall)) ///
					graphregion(color(white)) `gname' `twopt'
			}
		}
	}
	else {
		local ids ""
		frame `act' {
			gsort -expected
			forvalues j = 1/`=min(`topactors', `n')' {
				local ids "`ids' `=actorid[`j']'"
			}
		}
		local ntop : word count `ids'
		tempname pc
		frame copy `per' `pc'
		frame `pc' {
			quietly generate byte __in = 0
			foreach id of local ids {
				quietly replace __in = 1 if actorid == `id'
			}
			quietly collapse (sum) observed expected (min) tlo thi if __in, by(period)
			quietly generate double tmid = (tlo + thi) / 2
			format tlo thi %9.3f
			format observed %9.0f
			format expected %9.2f
			list period tlo thi observed expected, noobs sep(0) abbreviate(12)
			if "`noplot'" == "" {
				local gname = cond("`name'" != "", "name(`name', replace)", "")
				twoway (line expected tmid, lcolor(cranberry) lwidth(medthick)) ///
					(connected observed tmid, lcolor(navy) mcolor(navy) msize(small) lwidth(medthick)), ///
					xtitle("Time (mid point of the period)") ytitle("Number of events as sender") ///
					xlabel(, grid glcolor(gs14)) ylabel(, angle(0) grid glcolor(gs14) format(%9.0f)) ///
					legend(order(2 "Observed" 1 "Expected under the model") rows(1) size(small)) ///
					title("The `ntop' most active actors, events per period", size(medsmall)) ///
					graphregion(color(white)) `gname' `twopt'
			}
		}
		capture frame drop `pc'
	}
	return scalar N = `tobs'
	return scalar N_actors = `n'
	return scalar corr = `corr'
	return scalar pearson = `chi'
	return local over "`over'"
	if "`frame'" != "" {
		if "`over'" == "actor" frame copy `act' `frame'
		else frame copy `per' `frame'
		di as txt "Table in frame " as res "`frame'" as txt cond("`over'" == "actor", ": actorid, actor, observed, expected.", ": actorid, actor, period, tlo, thi, observed, expected (every actor and period).")
		return local frame "`frame'"
	}
	capture frame drop `act'
	capture frame drop `per'
end
