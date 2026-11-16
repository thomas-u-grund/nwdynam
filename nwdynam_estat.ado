*! version 1.8.0  16nov2026 nwdynam_estat
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
* The change point methods of estat changepoints are a port of the R package 'changepoint' (Rebecca Killick, Kaylea Haynes, Idris Eckley, Paul Fearnhead, Jamie Lee; GPL), see NOTICE.md.
* nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

* Postestimation of nwdynam: estat ic, estat eventlik, estat eventprob, estat outliers, estat changepoints.
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
