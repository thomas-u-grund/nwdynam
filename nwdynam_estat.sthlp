{smcl}
{* *! version 1.9.1  18nov2026 author: Thomas Grund}{...}
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
{p2col :{cmd:estat choiceplot}}rank of the realised receiver by the fitted probability (choice sub-models){p_end}
{p2col :{cmd:estat rateplot}}observed and expected events per actor and per period (rate sub-model){p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}{cmd:estat ic}{p_end}

{p 8 16 2}{cmd:estat eventlik} [{cmd:,} {opt fra:me(name)} {opt replace}]{p_end}

{p 8 16 2}{cmd:estat eventprob} {cmd:,} {opt fra:me(name)} [{opt replace} {opt maxr:ows(#)}]{p_end}

{p 8 16 2}{cmd:estat outliers} [{cmd:,} {opt met:hod(hampel|iqr|top)} {opt par:ameter(#)} {opt win:dow(#)} {opt fra:me(name)} {opt replace} {opt plot} {opt na:me(name)}]{p_end}

{p 8 16 2}{cmd:estat changepoints} [{cmd:,} {opt mom:ent(mean|variance)} {opt met:hod(pelt|amoc|binseg)} {opt win:dow(#)} {opt fra:me(name)} {opt replace} {opt plot} {opt na:me(name)}]{p_end}

{p 8 16 2}{cmd:estat choiceplot} [{cmd:,} {opt top:k(numlist)} {opt maxran:k(#)} {opt fra:me(name)} {opt replace} {opt nopl:ot} {opt na:me(name)} {opt maxr:ows(#)}]{p_end}

{p 8 16 2}{cmd:estat rateplot} [{cmd:,} {opt ov:er(actor|time)} {opt per:iods(#)} {opt topa:ctors(#)} {opt lab:el(#)} {opt fra:me(name)} {opt replace} {opt nopl:ot} {opt na:me(name)} {opt maxr:ows(#)}]{p_end}


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

{phang}{opt topk(numlist)} (estat choiceplot) lists the share of events whose realised receiver is among the k most probable candidates, for each k of the list (default 1 3 5 10).
{opt maxrank(#)} is the largest k of the graph (default 20; the curve ends where every risk set is covered). {opt noplot} gives the numbers without the graph.{p_end}

{phang}{opt over(actor|time)} (estat rateplot) chooses the graph: {bf:actor} (default) plots the observed against the expected number of events of every actor as sender with the 45 degree line, and labels
the {opt label(#)} (default 5) actors that deviate most; {bf:time} cuts the log into {opt periods(#)} (default 10) periods of equal length in time and plots the observed and expected number of events per
period of the {opt topactors(#)} (default 5) actors with the largest expected count.{p_end}

{phang}With {opt frame(name)} the commands store the per-event table (estat choiceplot) or the per-actor table, or the actor by period table with {cmd:over(time)} (estat rateplot). {opt maxrows()} is as in {cmd:estat eventprob}.{p_end}


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
{bf:Result plots (extensions).} goldfish has no such plots; the two commands are additions of this package and use the probabilities of {cmd:estat eventprob}, evaluated by the same native code.
{cmd:estat choiceplot} ranks, in every event, the realised receiver (for {bf:choice_coordination} the realised pair) among the members of the risk set by the fitted probability (rank 1 = most probable;
members with an equal probability share the mid rank, and a tie group that straddles k counts with the fraction of its places up to k). The graph shows the share of events with rank at most k, against
what a model without information would give (every member of the risk set equally probable) and against the expectation under the fitted model (the mean over the events of the sum of the k largest
probabilities); a fitted model whose observed curve lies close to the model curve and well above the uniform one describes who is chosen well. The output also gives the mean and the median rank.
{cmd:estat rateplot} adds up, over the events, the probability that each actor is the sender, giving the expected number of events of the actor; for the rate model with {opt intercept} the probability is the
actor's share of the total hazard, exp(lnrate) divided by the sum over the actors, so the expected counts always add up to the observed number of events (the plot checks who acts, not how many events occur).
Both commands need the sub-model they are made for (error 198 otherwise) and the fit in memory (error 498 as for the other commands).

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

{pstd}{cmd:estat choiceplot}: {bf:r(N)}, {bf:r(mean_rank)}, {bf:r(median_rank)}, {bf:r(mean_rank_uniform)}, {bf:r(mean_rank_relative)} (mean of rank divided by risk-set size), {bf:r(mean_riskset)}, the matrices {bf:r(topk)} (k, observed, uniform, model) and {bf:r(curve)} (the same for k = 1, 2, ...), {bf:r(frame)}. The frame has event, time, m (size of the risk set), rank and pobs.{p_end}

{pstd}{cmd:estat rateplot}: {bf:r(N)} (events), {bf:r(N_actors)}, {bf:r(corr)} (correlation of observed and expected counts), {bf:r(pearson)} (sum of (observed - expected)^2 / expected over the actors, descriptive only), {bf:r(over)}, {bf:r(frame)}. The frame has actorid, actor, observed and expected (with {cmd:over(time)} also period, tlo and thi, for every actor and period).{p_end}

{pstd}{cmd:estat changepoints}: {bf:r(n_cpts)}, {bf:r(N)}, {bf:r(window)}, {bf:r(method)}, {bf:r(moment)} and the row vector {bf:r(cpts)} with the last event of every regime (the last entry is always the number of events, as in the R package).{p_end}


{title:Also see}

{p 4 14 2}Help: {help nwdynam}{p_end}

{hline}
{pstd}{it:nwdynam version 1.9.1, 18 November 2026}{p_end}{* nwversion}
