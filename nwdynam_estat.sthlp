{smcl}
{* *! version 1.7.0  16nov2026 author: Thomas Grund}{...}
{marker topic}
{helpb nwtopical##analysis_dynamicmodels:[NW-2.6.8] Dynamic & Event-Based Models}

{title:Title}

{p2colset 9 29 31 2}{...}
{p2col :nwdynam postestimation {hline 2}}Postestimation tools for {help nwdynam}{p_end}
{p2colreset}{...}


{title:Description}

{pstd}
The following commands work after {cmd:nwdynam} for every sub-model. They offer what goldfish offers after estimation, with the same
definitions: the log likelihood of every event (goldfish option {cmd:returnIntervalLogL}), the probabilities of the events
(option {cmd:returnEventProbabilities}), {cmd:logLik}, {cmd:AIC} and {cmd:BIC}, {cmd:examineOutliers} and {cmd:examineChangepoints}.
They are evaluated at the coefficients in {cmd:e(b)} by the same native code as the fit. The data of the fit is kept in memory after
{cmd:nwdynam}; it is replaced by the next fit and gone after {cmd:clear all}, and the commands stop with return code 498 when the
estimation results on screen are not those of the fit in memory (for example after {cmd:estimates restore} of an earlier fit).

{p2colset 9 38 40 2}{...}
{p2col :{cmd:estat ic}}log likelihood, log likelihood per event, AIC and BIC{p_end}
{p2col :{cmd:estat eventlik}}log likelihood of every event{p_end}
{p2col :{cmd:estat eventprob}}model probabilities of every event{p_end}
{p2col :{cmd:estat outliers}}events with an unusually low log likelihood{p_end}
{p2col :{cmd:estat changepoints}}change points in the log likelihood of the event sequence{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}{cmd:estat ic}{p_end}

{p 8 16 2}{cmd:estat eventlik} [{cmd:,} {opt fra:me(name)} {opt replace}]{p_end}

{p 8 16 2}{cmd:estat eventprob} {cmd:,} {opt fra:me(name)} [{opt replace} {opt maxr:ows(#)}]{p_end}

{p 8 16 2}{cmd:estat outliers} [{cmd:,} {opt met:hod(hampel|iqr|top)} {opt par:ameter(#)} {opt win:dow(#)} {opt fra:me(name)} {opt replace} {opt plot} {opt na:me(name)}]{p_end}

{p 8 16 2}{cmd:estat changepoints} [{cmd:,} {opt mom:ent(mean|variance)} {opt met:hod(pelt|amoc|binseg)} {opt win:dow(#)} {opt fra:me(name)} {opt replace} {opt plot} {opt na:me(name)}]{p_end}


{title:Options}

{phang}{opt frame(name)} stores the table in a new frame; {opt replace} overwrites an existing frame. {cmd:estat eventprob} always needs a frame, because the
table has one row per event and candidate.{p_end}

{phang}{opt maxrows(#)} (estat eventprob) limits the number of rows of the table (events times candidates); the command stops with an error before it
builds a larger one. The default is 200000000.{p_end}

{phang}{opt method(), parameter(), window()} (estat outliers) choose the rule, as in goldfish: {bf:hampel} (the default; flags an event whose log likelihood
differs from the median of the {opt window()} events on each side by more than {opt parameter()} (default 3) scaled median absolute deviations; the default window is
half the number of events minus 1, as in goldfish, which tests only the middle event(s) of the log, so give {opt window()} for a rolling rule), {bf:iqr}
(flags log likelihoods below the median minus {opt parameter()}/2 interquartile ranges) and {bf:top} (flags the {opt parameter()} events with the lowest log likelihood).{p_end}

{phang}{opt moment(), method(), window()} (estat changepoints) choose what changes and how it is found: {bf:mean} (default) or {bf:variance} of the log likelihood; {bf:pelt}
(default), {bf:amoc} (at most one change) or {bf:binseg} (binary segmentation, up to 5 changes); {opt window()} is the minimum segment length in events (the default is the largest number of events that share a time point, as in
goldfish). The penalty is the modified Bayesian information criterion and the test statistic is the normal one.{p_end}

{phang}{opt plot} draws the log likelihood against time with the flagged events or the change points; {opt name()} names the graph.{p_end}


{title:Remarks}

{pstd}
{bf:Log likelihood of an event.} For {bf:choice} it is the log of the probability of the observed receiver; for {bf:rate} without intercept the log of the probability that the
observed sender acts; for {bf:choice_coordination} the log of the probability of the observed pair. In the rate model with intercept the log likelihood is a sum over the
periods between events: besides one record per event, goldfish includes the right-censored periods that a covariate change or the expiry of a window ends (they
are part of {cmd:e(ll)}). {cmd:estat eventlik} then reports one row per event and puts the censored periods before an event in the variable {cmd:llcens};
the part after the last event is in {cmd:r(ll_trailing)}. {cmd:estat ic} uses the number of records as goldfish does ({cmd:r(N)}), which is larger than the number of events in this case.
{cmd:estat outliers} and {cmd:estat changepoints} use the events only.

{pstd}
{bf:Probabilities.} Choice: one row per event and candidate receiver ({cmd:p}, summing to 1 over the candidates of an event); with composition change or {opt opportunities()} the candidates
are those present or available. Rate without intercept: one row per event and present actor, the probability that the actor is the sender. Rate with intercept: goldfish returns
for each actor the linear predictor minus the exposure of the period that ends at the event; this is not a probability and is stored as {cmd:lnrate}. Choice_coordination: one row per event
and pair of actors (actor1 < actor2); {cmd:p} is goldfish's matrix entry, twice the probability of the pair, so the entries of an event sum to 1.

{pstd}
{bf:Ties.} With {cmd:top}, events whose log likelihoods agree to the last digits are ranked in an arbitrary order; goldfish may then name another of them.

{pstd}
{bf:Credit.} The change point methods are a port of the R package {it:changepoint} (Rebecca Killick, Kaylea Haynes, Idris Eckley, Paul Fearnhead, Jamie Lee; GPL), which goldfish's
{cmd:examineChangepoints} calls; {bf:estat ic} and the other tools follow goldfish (Hollway, Stadtfeld, Hoffman, Uzaheta and contributors; GPL >= 3).


{title:Stored results}

{pstd}{cmd:estat ic} stores in {cmd:r()}: {bf:r(ll)}, {bf:r(ll_per_event)}, {bf:r(N)} (the records of the log likelihood, as goldfish counts), {bf:r(N_events)} (the events), {bf:r(df)}, {bf:r(aic)}, {bf:r(bic)} and the matrix {bf:r(S)}.{p_end}

{pstd}{cmd:estat eventlik}: {bf:r(N)}, {bf:r(ll_events)} (sum over the events), {bf:r(ll)}, {bf:r(ll_mean)}, {bf:r(ll_min)}, {bf:r(ll_max)}, {bf:r(n_cens)} and {bf:r(ll_trailing)} (rate model with intercept), {bf:r(frame)}.
The frame has event, origevent (the row of the event in the data), time, sender, receiver, ll (and llcens).{p_end}

{pstd}{cmd:estat eventprob}: {bf:r(N_rows)}, {bf:r(frame)}. The frame has event, time, actor and actorid (actor1, actor2 and their ids for choice_coordination), observed and {bf:p} ({bf:lnrate} for the rate model with intercept).{p_end}

{pstd}{cmd:estat outliers}: {bf:r(n_outliers)}, {bf:r(N)}, {bf:r(parameter)}, {bf:r(window)}, {bf:r(method)} and the row vector {bf:r(outliers)} of the event numbers in increasing order.{p_end}

{pstd}{cmd:estat changepoints}: {bf:r(n_cpts)}, {bf:r(N)}, {bf:r(window)}, {bf:r(method)}, {bf:r(moment)} and the row vector {bf:r(cpts)} with the last event of every regime (the last entry is always the number of events, as in the R package).{p_end}


{title:Also see}

{p 4 14 2}Help: {help nwdynam}{p_end}

{hline}
{pstd}{it:nwdynam version 1.7.0, 16 November 2026}{p_end}{* nwversion}
