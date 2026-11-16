/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
 * nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	ref_dynam_mata.do -- TEST-ONLY reference evaluators of nwdynam (Mata).

	These are the likelihood evaluators that nwdynam used before 1.5.0 (verified against goldfish). The package no longer fits with
	them: every fit runs in the plugin native/dynam_sim.c. The test suite loads this file after unw_dynam.do and compares the plugin
	with them (log likelihood, gradient, Hessian) on random small cases. Do not add this file to the library.
*/

capture mata: mata drop dynam_choice_eval_unit1()
capture mata: mata drop dynam_rate_eval_unit1()
capture mata: mata drop dynam_group_end()
capture mata: mata drop dynam_window_expire()
capture mata: mata drop dynam_choice_loglik_grad_unit1()
capture mata: mata drop dynam_choice_loglik_unit1()
capture mata: mata drop dynam_rate_loglik_grad_unit1()
capture mata: mata drop dynam_rate_loglik_unit1()
capture mata: mata drop dynam_choice_loglik_grad_multi()
capture mata: mata drop dynam_choice_loglik_multi()
capture mata: mata drop dynam_rate_loglik_grad_multi()
capture mata: mata drop dynam_rate_loglik_multi()
capture mata: mata drop dynam_rateint_loglik_grad_multi()
capture mata: mata drop dynam_rateint_loglik_multi()
capture mata: mata drop dynam_rateint_build_cache()
capture mata: mata drop dynam_rateint_cached_eval()
capture mata: mata drop dynam_coord_loglik_grad_multi()
capture mata: mata drop dynam_coord_loglik_multi()

mata:
mata set matastrict off


/*
	Simultaneous events (identical time stamps), as in goldfish: all events of a group are evaluated at the state
	as it was BEFORE the group; the state is updated only after the last event of the group. Returns the first
	row of the group if event i is the last of its group (apply the updates of rows gfirst..i now), else 0.
*/
real scalar dynam_group_end(pointer(class DynamState scalar) scalar pS, real scalar i) {
	if (i == (*pS).nevents) return((*pS).gfirst[i])
	if ((*pS).gfirst[i+1] == i+1) return((*pS).gfirst[i])
	return(0)
}

/*
	Windowed tie values (goldfish's window argument): W holds, for every pair, the summed value of the events already
	applied (rows 1..upto, chronological) whose time lies within `window' of tcur; ptr is the first such row. goldfish
	removes the increment of an event again at its time + window, after the dependent events of that time point, so an
	event exactly `window' before tcur still counts. sym = 1 updates both cells (undirected network).
*/
void dynam_window_expire(pointer(class DynamState scalar) scalar pS, real matrix W, real scalar ptr,
		real scalar upto, real scalar tcur, real scalar window, real scalar sym) {
	real scalar a, b
	while (ptr <= upto) {
		if (tcur - (*pS).events[ptr,3] <= window) break
		a = (*pS).events[ptr,1]
		b = (*pS).events[ptr,2]
		W[a,b] = W[a,b] - (*pS).events[ptr,4]
		if (sym) W[b,a] = W[b,a] - (*pS).events[ptr,4]
		ptr = ptr + 1
	}
}

/*
	DyNAM choice sub-model, v1 Unit 1: inertia + recip + indeg(alter),
	all unweighted/binary - three effect definitions read DIRECTLY from
	goldfish's own doc/goldfishEffects.Rmd (not guessed, not copied -
	see dev/dynam_unit1_crosscheck.R for the exact passages and the
	live R fit this was checked against):

	  inertia(s,j) = I(tie[s,j] > 0)   -- has s EVER sent to j before now?
	  recip(s,j)   = I(tie[j,s] > 0)   -- has j EVER sent to s before now?
	  indeg(j)     = colsum(tie[.,j])  -- j's own in-degree (count of
	                                       DISTINCT actors who have sent
	                                       to j before now) - "alter
	                                       type", per goldfishEffects.Rmd's
	                                       own type=c("alter","ego")
	                                       argument: the choice
	                                       sub-model always uses alter
	                                       type (candidate receiver's
	                                       own in-degree), never ego
	                                       type (that is the rate
	                                       sub-model's own use of
	                                       indeg(), not yet implemented
	                                       here)

	where tie[.,.] is a binary n x n matrix - NOT a raw event-count
	accumulator like unw_rem.do's own cij (a real, deliberate
	difference: goldfish's own default weighted=FALSE means REPEAT
	contacts do not inflate inertia/recip/indeg further, only the FIRST
	occurrence of a tie flips the corresponding cell from 0 to 1) -
	updated incrementally, STRICTLY prior to the current event (no
	lookahead), same "update AFTER using this event's own pre-event
	state" discipline as unw_rem.do's own cij/lastcontact.

	No first-event special-casing is needed (unlike unw_rem.do's own
	1/(n-1) fallback for its own degree-FRACTION effects, which divide
	by the number of prior events): tie starts all-zero, so
	inertia/recip/indeg are all simply 0 for every candidate at the
	first event - a well-defined value requiring no fallback branch,
	since these are counts/indicators, not fractions with an undefined
	0/0 at i=1.

	Risk set for event i: sender s_i is taken directly from the DATA
	(not modeled - that is the separate, not-yet-built rate sub-model,
	see this file's own header comment and the development notes).
	Every other actor j != s_i is a candidate receiver - a standard
	conditional (multinomial) logit over n-1 candidates, structurally
	the same softmax-over-risk-set machinery as unw_rem.do's own
	rem_loglik_grad_unit1(), but the risk set here is exactly ONE row
	(s_i's own row of a full n x n log-rate matrix), not the full n x n
	dyad space at once, since the sender is conditioned on rather than
	competing for the event itself.
*/
void dynam_choice_loglik_grad_unit1(real rowvector theta, pointer(class DynamState scalar) scalar pS,
		real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, mx, lrsum, g0, k
	real matrix tie, tieval
	real rowvector inertia_i, recip_i, indeg_i, lp, P

	n = (*pS).n
	ll = 0
	grad = J(1, 3, 0)
	tie = J(n, n, 0)
	tieval = J(n, n, 0)   // summed event values; a tie exists while the sum is above 0 (goldfish)

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]

		inertia_i = tie[s,.]
		recip_i = tie[.,s]'
		indeg_i = colsum(sign(tieval))   // goldfish: sign of each tie value (-1 for a value below 0)

		lp = theta[1] :* inertia_i :+ theta[2] :* recip_i :+ theta[3] :* indeg_i
		lp[s] = -1e300   // sender cannot choose itself - exclude from the risk set

		mx = max(lp)
		P = exp(lp :- mx)
		lrsum = mx + ln(sum(P))
		ll = ll + (lp[r] - lrsum)

		P = P :/ sum(P)   // proper softmax probability row over the risk set
		grad[1] = grad[1] + (inertia_i[r] - sum(inertia_i :* P))
		grad[2] = grad[2] + (recip_i[r] - sum(recip_i :* P))
		grad[3] = grad[3] + (indeg_i[r] - sum(indeg_i :* P))

		// update AFTER using this event's own pre-event state (after the whole group of equal times)
		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			tieval[(*pS).events[k,1], (*pS).events[k,2]] = tieval[(*pS).events[k,1], (*pS).events[k,2]] + (*pS).events[k,4]
			tie[(*pS).events[k,1], (*pS).events[k,2]] = (tieval[(*pS).events[k,1], (*pS).events[k,2]] > 0)
		}
	}
}

/*
	Pure log-likelihood wrapper (no gradient) - same role as
	unw_rem.do's own rem_loglik_unit1(), used by certification
	tests/diagnostics that only need the ll value at a given theta.
*/
real scalar dynam_choice_loglik_unit1(real rowvector theta, pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_choice_loglik_grad_unit1(theta, pS, ll, grad)
	return(ll)
}

/*
	DyNAM rate sub-model, Unit 2, NO-intercept case: indeg + outdeg,
	both "ego" type (an actor's OWN in/out degree - count of DISTINCT
	partners, unweighted - see this file's own header CORRECTION
	comment for how this formula was found, and dynam_choice_loglik_grad_unit1()'s
	own comment for why binary tie-presence rather than a raw event
	count is the right structure, same reasoning applies here).

	Risk set for event k: ALL n actors compete to be the next sender -
	unlike the choice sub-model (whose risk set is one already-known
	sender's own row, n-1 candidates), here every actor including the
	eventual non-senders is a candidate, so no "exclude self" step is
	needed (there is no "self" to exclude - the realized sender IS one
	of the n candidates, not compared against itself).

	  indeg_ego(i)  = colsum(tie[.,i])   -- count of distinct actors who have sent TO i before now
	  outdeg_ego(i) = rowsum(tie[i,.])   -- count of distinct actors i has sent TO before now

	ll_k = lp[s_k] - log(sum_i exp(lp[i]))   where lp[i] = theta[1]*indeg_ego(i) + theta[2]*outdeg_ego(i)

	Same incremental binary tie-matrix update as the choice sub-model
	("update AFTER using this event's own pre-event state") - in fact
	the exact same tie matrix, if both sub-models were ever fit
	together in one pass; kept as two independent functions here
	(matching goldfish's own separately-callable subModel="rate"/
	subModel="choice", and this file's own "two conditionally
	independent sub-models" framing) rather than sharing state across
	calls, since a caller may want either one alone.
*/
void dynam_rate_loglik_grad_unit1(real rowvector theta, pointer(class DynamState scalar) scalar pS,
		real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, mx, lrsum, g0, k
	real matrix tie, tieval
	real rowvector indeg_i, outdeg_i, lp, P

	n = (*pS).n
	ll = 0
	grad = J(1, 2, 0)
	tie = J(n, n, 0)
	tieval = J(n, n, 0)   // summed event values; a tie exists while the sum is above 0 (goldfish)

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]

		indeg_i = colsum(sign(tieval))   // goldfish: sign of each tie value (-1 for a value below 0)
		outdeg_i = rowsum(sign(tieval))'

		lp = theta[1] :* indeg_i :+ theta[2] :* outdeg_i

		mx = max(lp)
		P = exp(lp :- mx)
		lrsum = mx + ln(sum(P))
		ll = ll + (lp[s] - lrsum)

		P = P :/ sum(P)   // proper softmax probability row over the risk set (all n actors)
		grad[1] = grad[1] + (indeg_i[s] - sum(indeg_i :* P))
		grad[2] = grad[2] + (outdeg_i[s] - sum(outdeg_i :* P))

		// update AFTER using this event's own pre-event state (after the whole group of equal times)
		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			tieval[(*pS).events[k,1], (*pS).events[k,2]] = tieval[(*pS).events[k,1], (*pS).events[k,2]] + (*pS).events[k,4]
			tie[(*pS).events[k,1], (*pS).events[k,2]] = (tieval[(*pS).events[k,1], (*pS).events[k,2]] > 0)
		}
	}
}

real scalar dynam_rate_loglik_unit1(real rowvector theta, pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_rate_loglik_grad_unit1(theta, pS, ll, grad)
	return(ll)
}

/*
	Choice sub-model, generalized/selectable engine (the development notes'
	"effect selection" scope item - see this file's own header comment
	for how this relates to dynam_choice_loglik_grad_unit1() above).
	`active` is a 6-element 0/1 rowvector, fixed order (1=inertia,
	2=recip, 3=indeg, 4=same, 5=diff, 6=sim - the first three matching
	dynam_choice_loglik_grad_unit1()'s own effect order exactly; effects
	4-6 added 2026-09-02, "effect EXPANSION" - see this function's own
	header note below) - theta has sum(active) elements, in that same
	relative order. Effects 1-3: otherwise identical math to the unit-1
	function (same tie-matrix bookkeeping, same colsum(tie) computation
	of indeg_i - NOT the native port's own incremental-count
	optimization, since this is the MATA REFERENCE and should stay in
	the same form already verified against goldfish, not silently
	re-derive a different-but-equivalent formulation here too).

	Effects 4-6 (same/diff/sim - the development notes' "effect
	expansion" scope item, formulas read directly from goldfish's own
	doc/goldfishEffects.Rmd, not guessed): static per-actor ATTRIBUTE
	effects, not structural ones - a candidate's own value never changes
	across the fit, so (unlike inertia/recip/indeg) these need no
	incremental per-event state at all, just the sender's OWN attribute
	value broadcast against every candidate's:
	  same(s,j) = I(z_s = z_j)     -- homophily (same value)
	  diff(s,j) = |z_s - z_j|      -- heterophily (absolute difference)
	  sim(s,j)  = -|z_s - z_j|     -- homophily (similar value) - EXACTLY
	              -diff, by construction (goldfish's own documented
	              relationship, confirmed empirically: sim's coefficient
	              is diff's own coefficient with the sign flipped, same
	              |value|, same SE, same logLik - see
	              dev/dynam_unit4_attreffects_crosscheck.do). REAL
	              CONSEQUENCE, found while writing cscripts/test_nwdynam_ado.do:
	              diff() and sim() on the SAME variable are perfectly
	              collinear (one is an exact linear function of the
	              other), so a joint model including both is genuinely
	              UNIDENTIFIED - only the difference of their two
	              coefficients is estimable, not each separately. This
	              surfaced as a real optimizer failure ("flat region
	              encountered" / "simplex delta required"), not a
	              convergence-tuning problem to engineer around - using
	              two DIFFERENT variables for diff() and sim() together
	              is fine (they are then genuinely different dyadic
	              quantities); nwdynam.ado does not detect or reject the
	              same-variable case automatically (matching this
	              package's own general stance of not second-guessing
	              genuinely valid model specifications - a user might
	              deliberately want this contrast tested and hit the
	              expected non-identification directly).
	`same`/`diff`/`sim` each read from their OWN independent covariate
	vector (samevec/diffvec/simvec - NOT forced to share one variable,
	matching unw_rem.do's own covsndvec/covrecvec/covintvec convention
	of independent per-effect covariates - a real model may want a
	different variable for each). Static only, matching unw_rem.do's own
	v1 covariate scope (unit 4a). goldfish's own effect table confirms
	these apply to the CHOICE (and choice_coordination) sub-model only,
	never rate - no rate-side counterpart exists, matching
	dynam_rate_loglik_grad_multi()'s own unchanged 2-effect scope.

	Effect 7 (alter - "expansion batch 2", added 2026-09-02): the
	simplest possible attribute effect, per goldfish's own
	doc/goldfishEffects.Rmd - `s(i,j,t,z) = z_j`, the candidate
	receiver's OWN static covariate value, broadcast directly with no
	comparison to the sender at all (unlike same/diff/sim, which are
	genuinely dyadic). goldfish's own table confirms `alter` is
	choice-sub-model-only (never rate) - its rate-side counterpart is
	`ego` (effect 3 of dynam_rate_loglik_grad_multi() below, `r(i,t,z) =
	z_i`), a DIFFERENT effect name for what is really the same
	"static per-actor level" idea applied to the opposite sub-model's
	own risk-set role (candidate receiver vs. candidate sender) -
	matching goldfish's own naming exactly, not invented here.

	`windowinertia`/`windowrecip` ("effect expansion batch 3, added
	2026-09-02, user request "implement window effects," REDESIGNED
	same day as a per-effect suboption after user follow-up "check our
	goldfish syntax... make window a suboption of the effects... control
	this per effect" - `nwdynam.ado`'s own comment has the full account
	of why a single shared window scalar was replaced with two
	independent ones): a real-time recency cutoff on `inertia`/`recip`
	INDEPENDENTLY (goldfish's own `window=` argument is per-effect, not
	shared across a whole model - `inertia(net, window="7 days")` and
	`recip(net, window="1 day")` can genuinely differ in the same
	goldfish formula, and this engine now matches that exactly, not an
	earlier single-shared-window approximation of it). `indeg`/`outdeg`/
	`ego`/`alter` remain UNWINDOWED in this batch, a real, disclosed
	scope limit, not silently extended to match. >=1e300 means "no
	window" for that specific effect (the original always-sticky
	behavior, bit-identical to before this batch - see below); any
	finite value activates real-time filtering for THAT effect only.
	Formula, derived empirically and verified against real goldfish (not
	guessed - see dev/dynam_unit6_window_crosscheck.do and its own
	header comment for the toy-example-first derivation this project's
	own established discipline requires): a tie i->j counts as PRESENT
	at the current event's own evaluation time if and only if the MOST
	RECENT prior i->j event happened within that effect's own window of
	the current event's own timestamp (inclusive: `t_current -
	t_lastcontact <= window`) - NOT a decaying weight, a hard cutoff; a
	tie can EXPIRE (revert to absent) once its own most recent contact
	falls outside the window, even though it was present a moment
	before. This is a REAL DEPARTURE from every other effect in this
	file: the ordinal partial likelihood itself still uses only event
	ORDER (this file's own CORRECTION comment on the rate sub-model's
	own no-intercept case still holds - no waiting-time HAZARD term
	enters the likelihood), but window-eligible STATISTICS now read real
	elapsed time as a FILTER on which historical events still count,
	which the engine previously never touched at all -
	`(*pS).events[.,3]` (the time column DynamState already carries but
	effects 1-7 above never read) is used here for the first time.

	Implementation: the previous binary `tie` matrix is replaced by
	`lastcontact` (real time of the most recent i->j event, sentinel
	-1e300 for "never" - matching this file's own -1e300
	risk-set-exclusion sentinel convention). When both windows are
	inactive, `lastcontact[s,.] :> -1e300` reproduces the OLD `tie[s,.]`
	binary presence check exactly (bit-identical - verified directly,
	not merely argued, via dev/dynam_unit6_window_crosscheck.do's own
	regression check against the pre-window numbers); `indeg` (still
	unwindowed) is likewise derived from `lastcontact :> -1e300` rather
	than a separately-maintained matrix, so only ONE tracking structure
	is kept, not three redundant ones.

	UPDATE (event values): `lastcontact'/`tiecount' are replaced by `tieval', the value of each tie = the sum of the
	values of its events (column 4 of the event list, 1 without values), exactly as goldfish adds the increment of each
	event to the network. A tie exists while its value is above 0 (unweighted effects); weighted effects read the value.
	The unweighted indeg and outdeg add the SIGN of each tie value, as goldfish does: a tie whose value is below 0
	counts -1. Inside a window this happens when an event with a negative value outlives the event it cancelled.
	A window keeps its own matrix of the values of the events inside the window (dynam_window_expire(), goldfish's
	dissolve events at time + window), so `weighted' and a window combine as in goldfish. With every value positive all
	unweighted statistics are the same as before.

	Effect 8 (outdeg, "expansion batch 5", added 2026-09-02, user
	report that the help file's own "outdeg: rate sub-model only" label
	was checked against goldfish's own doc/goldfishEffects.Rmd and found
	WRONG - goldfish's own "Structural Effects" table lists outdeg as
	valid for BOTH rate and choice, not rate-only; this was a genuine
	missing effect, not a correct scope restriction): the choice
	sub-model's own "alter type" outdegree, `s(i,j,t,x) = sum_l
	I(x(t)_jl>0)` per goldfish's own doc - candidate j's own out-degree
	(how many distinct receivers j has themselves contacted), broadcast
	with no comparison to the sender, mirroring `indeg`'s own "alter
	type" shape exactly but transposed (`indeg` sums lastcontact's own
	COLUMN j - who has contacted j; `outdeg` sums its own ROW j - who j
	has contacted). Reuses the SAME `lastcontact` matrix already
	maintained for `inertia`/`recip`/`indeg` - no new tracking structure
	needed. Verified against real goldfish - see
	dev/dynam_unit8_choiceoutdeg_crosscheck.do. UNWINDOWED in this batch
	(matching `indeg`'s own still-unwindowed choice-side status,
	disclosed above) - a real, disclosed scope limit for a future batch,
	not silently extended to match `outdegwindow()`'s own existing
	RATE-side-only scope.

	The user separately asked, in the same exchange, whether
	`inertia`/`recip`/`same`/`diff`/`sim` (the OTHER effects the help
	file marks choice-only) could likewise be added to the rate
	sub-model. They cannot, and this is a structural fact about
	goldfish's own two sub-models, not an arbitrary scope choice:
	`inertia`/`recip`/`same`/`diff`/`sim` are all genuinely DYADIC
	statistics `s(i,j,t,x)` comparing a specific SENDER i against a
	specific CANDIDATE RECEIVER j (e.g. `recip`: does j have a recent
	tie back to i; `same`: do i and j share a covariate level) - the
	rate sub-model's own `r(i,t,x)` has no receiver role at all (it is a
	single-actor hazard over who acts next), so a dyadic comparison has
	no argument to compare i against. `indeg`/`outdeg` are the exception
	precisely because they are NOT dyadic - each is a property of a
	SINGLE actor (in-degree or out-degree), expressible identically
	whether that actor is playing the rate sub-model's own "candidate
	sender" role or the choice sub-model's own "candidate receiver"
	role. `alter`'s own single-actor analogue already exists under
	goldfish's own different name for the rate sub-model - `ego` (effect
	3 of dynam_rate_loglik_grad_multi() below) - so there is no missing
	"alter for rate" to add either; it was already implemented under its
	goldfish-given name in expansion batch 2.

	Effects 9-14 ("expansion batch 6", added 2026-09-02, user: "work
	your way through the list" - a systematic audit against goldfish's
	FULL effect catalog, not just the ones already covered): the
	two-path closure effects (`trans`, `cycle`, `commonSender`,
	`commonReceiver`, `four`) and `nodeTrans`, per goldfish's own
	doc/goldfishEffects.Rmd formulas, computed on the DEPENDENT network
	itself (goldfish's own stated default when no explicit network
	argument is given - "For network effects, this would be either an
	explanatory network or the dependent network (the default)"). Using
	an EXOGENOUS network instead (goldfish's own `trans(othernet)`, or
	the two-network `mixedTrans`/`mixedCycle`/`mixedCommonSender`/
	`mixedCommonReceiver` variants) is a real, disclosed scope limit for
	a future batch - nwdynam has no mechanism yet for an effect to read
	a network other than the one being predicted.

	All five closure effects are two-path counts over the SAME binary
	presence matrix `A` (`lastcontact :> -1e300`, reusing the existing
	tracking structure exactly like `indeg`/`outdeg` already do) -
	computed via matrix products rather than a per-candidate loop, since
	each is a full ROW (or column) of a single n x n product:
	  trans_i(j)          = sum_k A[i,k]A[k,j]        = (A*A)[i,.]
	  cycle_i(j)           = sum_k A[j,k]A[k,i]        = (A*A)[.,i]'
	  commonSender_i(j)    = sum_k A[k,i]A[k,j]        = (A'*A)[i,.]
	  commonReceiver_i(j)  = sum_k A[i,k]A[j,k]        = (A*A')[i,.]
	  four_i(j)            = sum_kl A[i,k]A[l,k]A[l,j], k,l NOT IN {i,j}
	                           and k != l (see correction below)
	with `i` bound to the realized sender `s` at each event (each is
	genuinely dyadic - varies by candidate `j`, matching `inertia`'s own
	sender-dependent shape, not `indeg`'s own candidate-only shape).

	`four` needed a REAL correction beyond the naive matrix product,
	found by direct comparison against goldfish's own reconstructed
	statistic matrix (via `estimate(..., preprocessingOnly=TRUE)`'s own
	`dependentStatsChange` deltas), not assumed from the documented
	formula alone: the naive `(A*(A'*A))[i,.]` (matching trans/cycle/
	commonSender/commonReceiver's own pattern exactly) reproduced
	trans/cycle/commonSender/commonReceiver bit-for-bit but was
	systematically almost DOUBLE goldfish's own real `four` values and
	produced a wildly wrong fitted coefficient (1.644 vs goldfish's own
	0.897, logLik off by 660) - `four`'s own doc/goldfishEffects.Rmd
	prose ("closes more three-paths i->k<-l->j") requires i, j, k, l to
	be four GENUINELY DISTINCT actors, unlike the other four closure
	effects' own two-path formulas (which only involve 3 roles, i/j/k,
	where the existing zero-diagonal convention already excludes every
	degenerate case automatically). The closed-form correction, derived
	and verified via a small hand-traceable random toy network (n=6)
	against a full four-nested-loop brute-force distinct-actors
	computation BEFORE touching this file, then confirmed on real data
	(matching goldfish's own reconstructed statistic matrix off-diagonal
	to within rounding):
	  four_i(j) = (A*(A'*A))[i,j] - A[i,j] * (rowsum(A[i,.]) +
	              diagonal(A'*A)[j] - A[i,j])
	subtracting exactly the k=i, k=j, and l=i contributions from the
	naive sum (l=j and k=k self-loop contributions are already zero via
	the existing zero-diagonal convention, so only these three terms
	needed an explicit correction). Re-verified against real goldfish
	after the fix - see dev/dynam_unit9_closure_crosscheck.do - matching
	to within 8e-6 on the coefficient, 0.02 on logLik.
	`nodeTrans` (the one BOTH-submodel effect in this batch, per
	goldfish's own table) is different in kind: `r(i) = sum_j A[i,j] *
	(A*A)[i,j]`, a property of a SINGLE actor with no dependence on the
	realized sender at all (matching `indeg`/`outdeg`'s own
	candidate-only shape, computed once per event as `rowsum(A :* (A*A))`
	and reused directly for every candidate/actor).
	`A*A`/`A'*A`/`A*A'` are each computed CONDITIONALLY, only when an
	effect that needs them is active (`needAA`/`needAtA`/`needAAt`), to
	avoid the O(n^3) matrix-multiply cost when none of these 6 effects
	are requested - every other existing effect remains exactly as cheap
	as before this batch. UNWINDOWED in this batch, a real disclosed
	scope limit (goldfish's own table shows `window` IS valid on these
	effects; matching the already-disclosed unwindowed status of
	choice's own `indeg`, not silently extended to match). Verified
	against real goldfish on Social_Evolution - see
	dev/dynam_unit9_closure_crosscheck.R/.do, including a genuine
	multi-effect model (`inertia + recip + trans + cycle + commonSender
	+ commonReceiver` together) matching to within 0.01 on every
	coefficient. `commonReceiver` triggers a real goldfish-side
	`isTwoMode = TRUE` auto-detection WARNING (not an error) on this
	directed dataset - a goldfish-internal heuristic quirk, not caused
	by anything in nwdynam - checked directly and found NOT to change
	goldfish's own reported coefficient from what the directed,
	isTwoMode=FALSE formula above predicts (matched to within 0.001).

	Effects 15-16 ("expansion batch 7", added 2026-09-02, continuing
	"work your way through the list"): `egoAlterInt` and `tertius`, per
	goldfish's own doc/goldfishEffects.Rmd formulas.

	Effect 15 (egoAlterInt): `s(i,j,t,z1,z2) = z1_i * z2_j` - the
	SENDER's own covariate 1 value multiplied by the CANDIDATE's own
	covariate 2 value (an interaction, not a comparison - the two
	covariates can be the same variable or different ones, matching
	goldfish's own two-argument `egoAlterInt(attribute1, attribute2)`).
	Choice-only per goldfish's own table (needs both a sender role and a
	candidate-receiver role - the same structural reason `inertia`/
	`recip`/`same`/`diff`/`sim` are choice-only, see the earlier comment
	in this same header).

	Effect 16 (tertius, "alter" type for choice / "ego" type for rate -
	the SAME underlying computation reused across both engines exactly
	like `nodeTrans`, since it is a property of a single actor with no
	dependence on the realized sender): goldfish's own formula -
	`s(i,j,t,x,z) = mean_{k: x_kj>0}(z_k)` - the MEAN covariate value of
	the candidate's own in-neighbors (senders who have contacted them),
	on the DEPENDENT network itself (goldfish's own default network
	argument, matching the closure-effect batch's own scope decision).
	v1 scope, disclosed: only goldfish's own DEFAULT `aggregateFun`
	(mean) is supported - `sum`/`max`/other custom aggregation functions
	are not requestable. Computed via `tertiusnumer = tertiusvec * A`
	(matrix product: for each candidate j, sums the covariate values of
	everyone who has contacted them) and `tertiusdenom = indeg_i`
	(already computed each event, reused directly rather than a
	redundant `colsum(A)`), with an explicit `hasin' mask to avoid 0/0
	propagating as a Mata missing value.

	ISOLATE IMPUTATION - a real, documented-vs-actual discrepancy found
	by direct comparison, not trusted from the doc alone (same
	discipline as the `four()` fix above): goldfish's own doc/
	goldfishEffects.Rmd states "When a node does not have in-neighbors,
	the tertius effect is imputed as the average of the aggregate
	values of nodes with in-neighbors." Implementing that literally
	first (a `hasin'-weighted global mean of every non-isolated actor's
	own tertius value) produced a coefficient wildly different from
	real goldfish (1.153 off, logLik off by 151 on `choice ~ tertius`
	alone). `estimate(..., preprocessingOnly=TRUE)`'s own
	`dependentStatsChange` deltas were reconstructed into goldfish's own
	actual statistic matrix for BOTH sub-models (same technique as the
	`four()` investigation) and compared directly against the
	global-mean-imputed values: every non-isolated actor matched
	EXACTLY, but every isolated actor's own real value was 0, not the
	documented global mean - confirmed on both the choice (candidate
	role) and rate (single-actor role) sub-models independently, ruling
	out a submodel-specific explanation. The documented "Note" does not
	match goldfish 1.6.12's own actual runtime behavior for this
	dataset - isolates are imputed as 0, full stop, matching this
	engine's own pre-existing sentinel convention for "nothing here yet"
	elsewhere in the file. Implemented as verified, not as documented.
	Re-verified against real goldfish on Social_Evolution after the fix
	- see dev/dynam_unit10_tertius_egoalterint_crosscheck.R/.do -
	matching to within 0.01 on every coefficient, including `tertius`
	under BOTH sub-models and a genuine 3-effect combined choice model
	(`inertia + alter + egoAlterInt`).

	`n1' ("expansion batch 9", two-mode/bipartite support, added
	2026-09-02, user: "continue" - the user's own explicit go-ahead
	after being shown the remaining large architectural items).
	Previously rejected outright by nwdynam.ado - see that file's own
	updated comment for why the guard's stated reason no longer holds
	and what was verified before lifting it.

	MECHANISM: `n1' is the mode-1 actor count (0 = not two-mode). This
	package's own combined actor-index convention (established by
	nwset's own bipartite ingestion, already reused by nwergm's own
	bipartite() support) puts mode-1 actors at indices 1..n1 and mode-2
	at n1+1..n, contiguous. goldfish's own two-mode DyNAM architecture
	assumes a STRICTLY ONE-DIRECTIONAL bipartite dependent network -
	confirmed directly via a hand-built toy affiliation network (people
	sending to orgs, never the reverse) - every event's sender is drawn
	from node set 1 ("nodes"), every candidate/receiver from node set 2
	("nodes2"), never mixed. Implemented as a simple risk-set MASK on
	top of the existing engine: choice sub-model excludes mode-1 actors
	from the candidate set (`lp[1..n1] = -1e300`); rate sub-model
	excludes mode-2 actors from the "who acts next" risk set
	(`lp[(n1+1)..n] = -1e300`). Hand-traced against goldfish's own
	`dependentStatsChange` deltas (same `preprocessingOnly=TRUE'
	technique used for `four()'/`tertius()' above) for `inertia' on a
	tiny 2-person/2-org/4-event example FIRST, confirming the naive
	"just restrict the risk set, formula unchanged" hypothesis holds for
	genuinely dyadic tie-presence effects - only THEN extended to a
	larger (6-person/4-org/60-event) toy dataset for full verification.

	PER-EFFECT ELIGIBILITY - goldfish's own documented `isTwoMode`
	column (doc/goldfishEffects.Rmd) turned out UNRELIABLE as a
	predictor of which effects actually work for two-mode - checked
	directly, not trusted: several effects the table marks eligible
	(`recip`, `outdeg` choice-side/alter-type, `commonReceiver`, `indeg`
	rate-side/ego-type) are hard-REJECTED by goldfish's own engine at
	runtime, each for the SAME real structural reason - goldfish's own
	one-directional architecture means mode 1 has no in-ties and mode 2
	has no out-ties, so any effect needing that missing role is
	undefined, not just usually-zero. These rejections are MIRRORED
	here (not silently allowed to produce a degenerate fit).

	A DEEPER finding, going beyond simple risk-set restriction: a few
	effects goldfish DOES accept for two-mode (`commonSender`, `four`)
	produced NON-degenerate fits that pure risk-set-restriction
	reasoning could not explain, matching goldfish's own doc warning
	("we omit the difference in the computation of the statistics when
	isTwoMode is used") - meaning some effects genuinely redefine their
	own formula for two-mode, not just their risk set. This was
	confirmed to run deeper than expected: direct comparison against
	goldfish's own reconstructed statistics found `same()` and `ego()` -
	the SIMPLEST possible attribute effects, a plain value comparison
	and a plain value broadcast respectively - ALSO silently wrong under
	the naive "read one combined covariate vector spanning both modes"
	approach originally attempted here (`same()` off by a large margin;
	`ego()` off by a smaller but still real margin). Investigating why
	revealed goldfish's own internal two-mode statistic array applies
	its ordinary "zero the self-tie" convention by matching raw ROW
	INDEX to raw COLUMN INDEX, even though rows and columns span two
	GENUINELY DIFFERENT node sets in a two-mode network (so
	person-index-3 and org-index-3 get treated as "the same actor" for
	self-exclusion purposes, despite being unrelated actors in different
	modes) - a real goldfish-internal representation quirk that could
	not be confidently reverse-engineered into a correct general formula
	in the time available. Given this, EVERY attribute effect (`same`/
	`diff`/`sim`/`alter`/`ego`/`egoAlterInt`/`tertius`) is rejected for
	two-mode in nwdynam.ado, not just the ones goldfish itself hard-
	rejects - a disclosed, unresolved gap, not a silently-shipped guess.

	VERIFIED WORKING, with an EXACT numerical match to real goldfish on
	the 6-person/4-org/60-event toy dataset (dev/dynam_unit12_twomode_
	crosscheck.R/.do): `inertia` (choice: -0.42729, alone; -0.73342
	alongside `indeg`), `indeg` alter-type (choice: 0.12673 alone;
	0.23018 alongside `inertia`), `outdeg` ego-type (rate: 0.23797) -
	matching to within 1e-5 on every coefficient. `nodeTrans`/`trans`/
	`cycle` are excluded per goldfish's own `isTwoMode=×` table entries
	(not independently re-verified, since every OTHER table entry that
	WAS re-verified turned out unreliable in at least one direction -
	disclosed, not assumed safe either way). `commonSender`/`four`
	(goldfish accepts them, but with the unexplained formula difference
	above) and windowed/weighted effects (not yet verified together with
	two-mode at all) remain real, disclosed gaps for a future batch.
*/
void dynam_choice_loglik_grad_multi(real rowvector theta, real rowvector active,
		real matrix covmat, real rowvector windowvec, real rowvector weightvec, real scalar n1,
		real matrix oppmat, real matrix tiemat,
		pointer(class DynamState scalar) scalar pS, real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, mx, lrsum, j, e, tcur, nowindowinertia, nowindowrecip
	real scalar needAA, needAtA, needAAt, windowinertia, windowrecip
	real scalar weightedinertia, weightedrecip, weightedindeg, weightedoutdeg, g0, k, a, b, v, pwi, pwr
	real matrix tieval, winertia, wrecip, A, AA, AtA, AAt, AAtA
	real rowvector inertia_i, recip_i, indeg_i, outdeg_i, same_i, diff_i, sim_i
	real rowvector trans_i, cycle_i, commonsender_i, commonreceiver_i, four_i, nodetrans_i
	real rowvector egoalterint_i, tertius_i, tertiusnumer, tertiusdenom, hasin, tie_i
	real rowvector lp, P, gradfull, idx
	real rowvector samevec, diffvec, simvec, altervec, tertiusvec, egoalterint1vec, egoalterint2vec

	// `covmat' bundles every per-actor covariate vector into ONE 7 x n
	// matrix (row 1=same, 2=diff, 3=sim, 4=alter, 5=tertius,
	// 6=egoalterint1, 7=egoalterint2), `windowvec' bundles both window
	// scalars into ONE 1 x 2 vector (1=windowinertia, 2=windowrecip),
	// and `weightvec' bundles the four `weighted=TRUE' modifier flags
	// into ONE 1 x 4 vector (1=inertia, 2=recip, 3=indeg, 4=outdeg) -
	// rather than passing each as its own optimize_init_argument()
	// slot. Stata's own optimize() framework caps optimize_init_argument()
	// at index 9 (confirmed directly by a minimal repro, not assumed -
	// the previous 11-separate-argument design for this function hit
	// that cap exactly, silently (`invalid argument index'), once
	// egoAlterInt/tertius pushed the count past it in "expansion batch
	// 7"). Bundling into matrices instead of separate arguments keeps
	// this function's own outer signature FIXED regardless of how many
	// more per-actor covariates or modifier flags a future batch adds -
	// only `covmat'/`weightvec' grow, no new optimize_init_argument()
	// slot ever needed again.
	samevec = covmat[1,.]
	diffvec = covmat[2,.]
	simvec = covmat[3,.]
	altervec = covmat[4,.]
	tertiusvec = covmat[5,.]
	egoalterint1vec = covmat[6,.]
	egoalterint2vec = covmat[7,.]
	windowinertia = windowvec[1]
	windowrecip = windowvec[2]
	weightedinertia = weightvec[1]
	weightedrecip = weightvec[2]
	weightedindeg = weightvec[3]
	weightedoutdeg = weightvec[4]

	n = (*pS).n
	ll = 0
	gradfull = J(1, 17, 0)
	// tieval: the value of every tie, the sum of the values of its events so far (goldfish: each event adds its value,
	// the increment); a tie exists while its value is above 0 (A). winertia/wrecip: the same within the windows.
	tieval = J(n, n, 0)
	winertia = J(n, n, 0)
	wrecip = J(n, n, 0)
	pwi = 1
	pwr = 1
	same_i = J(1, n, 0)
	diff_i = J(1, n, 0)
	sim_i = J(1, n, 0)
	nowindowinertia = (windowinertia >= 1e300)
	nowindowrecip = (windowrecip >= 1e300)
	needAA = active[9] | active[10] | active[14]
	needAtA = active[11] | active[13]
	needAAt = active[12]
	if (active[15]) egoalterint_i = J(1, n, 0)

	idx = J(1, 17, 0)
	j = 0
	for (e=1; e<=17; e++) {
		if (active[e]) {
			j = j + 1
			idx[e] = j
		}
	}

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]
		tcur = (*pS).events[i,3]

		A = (tieval :> 0)
		if (!nowindowinertia) dynam_window_expire(pS, winertia, pwi, (*pS).gfirst[i] - 1, tcur, windowinertia, 0)
		if (!nowindowrecip) dynam_window_expire(pS, wrecip, pwr, (*pS).gfirst[i] - 1, tcur, windowrecip, 0)
		if (nowindowinertia) {
			if (weightedinertia) inertia_i = tieval[s,.]
			else inertia_i = A[s,.]
		}
		else {
			if (weightedinertia) inertia_i = winertia[s,.]
			else inertia_i = (winertia[s,.] :> 0)
		}
		if (nowindowrecip) {
			if (weightedrecip) recip_i = tieval[.,s]'
			else recip_i = A[.,s]'
		}
		else {
			if (weightedrecip) recip_i = wrecip[.,s]'
			else recip_i = (wrecip[.,s]' :> 0)
		}
		if (weightedindeg) indeg_i = colsum(tieval)
		else indeg_i = colsum(sign(tieval))
		if (active[8]) {
			if (weightedoutdeg) outdeg_i = rowsum(tieval)'
			else outdeg_i = rowsum(sign(tieval))'
		}
		if (active[4]) same_i = (samevec :== samevec[s])
		if (active[5]) diff_i = abs(diffvec :- diffvec[s])
		if (active[6]) sim_i = -abs(simvec :- simvec[s])
		if (active[15]) egoalterint_i = egoalterint1vec[s] :* egoalterint2vec

		if (needAA) AA = A * A
		if (needAtA) AtA = A' * A
		if (needAAt) AAt = A * A'
		if (active[9]) trans_i = AA[s,.]
		if (active[10]) cycle_i = AA[.,s]'
		if (active[11]) commonsender_i = AtA[s,.]
		if (active[12]) commonreceiver_i = AAt[s,.]
		if (active[13]) {
			AAtA = A * AtA
			four_i = AAtA[s,.] :- A[s,.] :* (rowsum(A[s,.]) :+ diagonal(AtA)' :- A[s,.])
		}
		if (active[14]) nodetrans_i = rowsum(A :* AA)'
		if (active[16]) {
			tertiusnumer = tertiusvec * A
			// Always the count of DISTINCT in-neighbors, never the
			// weighted indegree - `indeg_i' above may itself be
			// weighted (`weightedindeg') if `indeg' is ALSO active in
			// the SAME model, which must not silently change tertius's
			// own denominator.
			tertiusdenom = colsum(A)
			hasin = (tertiusdenom :> 0)
			tertius_i = hasin :* (tertiusnumer :/ (tertiusdenom :+ (1 :- hasin)))
		}

		// tie() ("cross-network effects, v1 scope", added 2026-09-02):
		// unweighted presence in a SEPARATE, STATIC exogenous network
		// (`tiemat', an n x n 0/1 matrix, unchanging across events - not
		// itself an event-declared/evolving network, a real, disclosed
		// v1 scope limit) - `s(i,j,t,x) = I(x_ij>0)' per goldfish's own
		// doc/goldfishEffects.Rmd, identical in shape to `inertia' but
		// reading a FIXED matrix instead of the incrementally-tracked
		// `lastcontact'. Verified directly against real goldfish (see
		// dev/dynam_unit18_tie_crosscheck.R/.do) - `tiemat' has 0 rows
		// when tie() was not given (matching `oppmat''s own "not given"
		// sentinel convention), skipping the row-extraction entirely.
		if (active[17]) tie_i = tiemat[s,.]

		lp = J(1, n, 0)
		if (active[1]) lp = lp :+ theta[idx[1]] :* inertia_i
		if (active[2]) lp = lp :+ theta[idx[2]] :* recip_i
		if (active[3]) lp = lp :+ theta[idx[3]] :* indeg_i
		if (active[4]) lp = lp :+ theta[idx[4]] :* same_i
		if (active[5]) lp = lp :+ theta[idx[5]] :* diff_i
		if (active[6]) lp = lp :+ theta[idx[6]] :* sim_i
		if (active[7]) lp = lp :+ theta[idx[7]] :* altervec
		if (active[8]) lp = lp :+ theta[idx[8]] :* outdeg_i
		if (active[9]) lp = lp :+ theta[idx[9]] :* trans_i
		if (active[10]) lp = lp :+ theta[idx[10]] :* cycle_i
		if (active[11]) lp = lp :+ theta[idx[11]] :* commonsender_i
		if (active[12]) lp = lp :+ theta[idx[12]] :* commonreceiver_i
		if (active[13]) lp = lp :+ theta[idx[13]] :* four_i
		if (active[14]) lp = lp :+ theta[idx[14]] :* nodetrans_i
		if (active[15]) lp = lp :+ theta[idx[15]] :* egoalterint_i
		if (active[16]) lp = lp :+ theta[idx[16]] :* tertius_i
		if (active[17]) lp = lp :+ theta[idx[17]] :* tie_i
		lp[s] = -1e300
		// Two-mode: mode-1 actors (indices 1..n1, which already
		// includes the sender itself, so this subsumes the line above)
		// can never be a CANDIDATE - goldfish's own architecture
		// assumes mode 2 only ever receives (see nwdynam.ado's own
		// two-mode comment for the real-goldfish verification this
		// mirrors).
		if (n1 > 0) lp[1..n1] = J(1, n1, -1e300)
		// opportunitiesList (goldfish's own estimationInit argument,
		// choice sub-model only): per-event candidate restriction -
		// actors NOT marked available for event i are excluded from the
		// risk set for THAT event only, matching goldfish's own
		// behavior exactly (verified directly: the sender's own
		// self-exclusion above still applies unconditionally regardless
		// of whether the sender itself is listed as "available").
		// `oppmat' has 0 rows when opportunitiesList was not given (the
		// sentinel for "no restriction," skipping this entirely).
		if (rows(oppmat) > 0) lp = lp :+ (oppmat[i,.] :== 0) :* (-1e300)

		mx = max(lp)
		P = exp(lp :- mx)
		lrsum = mx + ln(sum(P))
		ll = ll + (lp[r] - lrsum)

		P = P :/ sum(P)
		if (active[1]) gradfull[1] = gradfull[1] + (inertia_i[r] - sum(inertia_i :* P))
		if (active[2]) gradfull[2] = gradfull[2] + (recip_i[r] - sum(recip_i :* P))
		if (active[3]) gradfull[3] = gradfull[3] + (indeg_i[r] - sum(indeg_i :* P))
		if (active[4]) gradfull[4] = gradfull[4] + (same_i[r] - sum(same_i :* P))
		if (active[5]) gradfull[5] = gradfull[5] + (diff_i[r] - sum(diff_i :* P))
		if (active[6]) gradfull[6] = gradfull[6] + (sim_i[r] - sum(sim_i :* P))
		if (active[7]) gradfull[7] = gradfull[7] + (altervec[r] - sum(altervec :* P))
		if (active[8]) gradfull[8] = gradfull[8] + (outdeg_i[r] - sum(outdeg_i :* P))
		if (active[9]) gradfull[9] = gradfull[9] + (trans_i[r] - sum(trans_i :* P))
		if (active[10]) gradfull[10] = gradfull[10] + (cycle_i[r] - sum(cycle_i :* P))
		if (active[11]) gradfull[11] = gradfull[11] + (commonsender_i[r] - sum(commonsender_i :* P))
		if (active[12]) gradfull[12] = gradfull[12] + (commonreceiver_i[r] - sum(commonreceiver_i :* P))
		if (active[13]) gradfull[13] = gradfull[13] + (four_i[r] - sum(four_i :* P))
		if (active[14]) gradfull[14] = gradfull[14] + (nodetrans_i[r] - sum(nodetrans_i :* P))
		if (active[15]) gradfull[15] = gradfull[15] + (egoalterint_i[r] - sum(egoalterint_i :* P))
		if (active[16]) gradfull[16] = gradfull[16] + (tertius_i[r] - sum(tertius_i :* P))
		if (active[17]) gradfull[17] = gradfull[17] + (tie_i[r] - sum(tie_i :* P))

		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			a = (*pS).events[k,1]
			b = (*pS).events[k,2]
			v = (*pS).events[k,4]
			tieval[a,b] = tieval[a,b] + v
			if (!nowindowinertia) winertia[a,b] = winertia[a,b] + v
			if (!nowindowrecip) wrecip[a,b] = wrecip[a,b] + v
		}
	}

	grad = J(1, 0, 0)
	for (e=1; e<=17; e++) if (active[e]) grad = (grad, gradfull[e])
}

real scalar dynam_choice_loglik_multi(real rowvector theta, real rowvector active,
		real matrix covmat, real rowvector windowvec, real rowvector weightvec, real scalar n1,
		real matrix oppmat, real matrix tiemat,
		pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_choice_loglik_grad_multi(theta, active, covmat, windowvec, weightvec, n1, oppmat, tiemat, pS, ll, grad)
	return(ll)
}

/*
	Rate sub-model, generalized/selectable engine - same role as
	dynam_choice_loglik_grad_multi() above, for the rate sub-model's own
	effects (1=indeg, 2=outdeg, matching dynam_rate_loglik_grad_unit1()'s
	own order; 3=ego, added 2026-09-02 - "expansion batch 2").

	Effect 3 (ego): the simplest possible attribute effect, per
	goldfish's own doc/goldfishEffects.Rmd - `r(i,t,z) = z_i`, the
	candidate (potential next sender)'s OWN static covariate value,
	broadcast directly with no comparison at all - the rate-side
	counterpart of the choice sub-model's own `alter` effect
	(dynam_choice_loglik_grad_multi()'s effect 7), matching goldfish's
	own naming exactly (`ego` for rate, `alter` for choice - the SAME
	underlying "static per-actor level" idea, just named for whichever
	risk-set role that sub-model's own candidates play). goldfish's own
	effect table confirms `ego` is rate-sub-model-only (never choice).

	`windowindeg`/`windowoutdeg` (added 2026-09-02, user follow-up "can
	we make [window] also available for the rate part of the model?",
	immediately after the choice-side inertiawindow()/recipwindow()
	self-activation redesign): a real-time recency cutoff on indeg/outdeg
	INDEPENDENTLY, exactly mirroring windowinertia/windowrecip's own
	semantics and >=1e300 "no window" sentinel. goldfish's own
	doc/goldfishEffects.Rmd effect table lists `window` as a real
	argument on indeg/outdeg and confirms both apply to the RATE
	sub-model (not just choice). Formula (same shape as the choice
	engine's own inertia/recip window, re-derived and verified against
	real goldfish here rather than assumed identical - see
	dev/dynam_unit7_ratewindow_crosscheck.R/.do): windowed indeg for
	candidate i counts only senders j whose MOST RECENT prior j->i event
	happened within the window of the current event's own timestamp;
	windowed outdeg for candidate i counts only receivers j whose most
	recent prior i->j event is within the window - the same
	`t_current - t_lastcontact <= window` hard cutoff, not a decaying
	weight. Implementation mirrors the choice engine exactly: the
	previous binary `tie` matrix is replaced by `lastcontact` (real time
	of the most recent i->j event, sentinel -1e300 for "never"). When
	both windows are inactive, `lastcontact :> -1e300` reproduces the
	OLD `tie :> 0` binary presence check exactly (bit-identical -
	verified directly via dev/dynam_unit7_ratewindow_crosscheck.do's own
	no-window regression sub-check against unit1's pre-window numbers).
*/
void dynam_rate_loglik_grad_multi(real rowvector theta, real rowvector active,
		real rowvector egovec, real rowvector tertiusvec, real scalar windowindeg, real scalar windowoutdeg,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1,
		pointer(class DynamState scalar) scalar pS, real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, mx, lrsum, j, e, tcur, nowindowindeg, nowindowoutdeg, g0, k, a, b, v, pwd, pwo
	real matrix tieval, windeg, woutdeg, A, AA
	real rowvector indeg_i, outdeg_i, nodetrans_i, tertius_i, tertiusnumer, tertiusdenom, hasin, lp, P, gradfull, idx

	n = (*pS).n
	ll = 0
	gradfull = J(1, 5, 0)
	// tieval, windeg, woutdeg: tie values (summed event values) overall and within the windows - see
	// dynam_choice_loglik_grad_multi()
	tieval = J(n, n, 0)
	windeg = J(n, n, 0)
	woutdeg = J(n, n, 0)
	pwd = 1
	pwo = 1
	nowindowindeg = (windowindeg >= 1e300)
	nowindowoutdeg = (windowoutdeg >= 1e300)

	idx = J(1, 5, 0)
	j = 0
	for (e=1; e<=5; e++) {
		if (active[e]) {
			j = j + 1
			idx[e] = j
		}
	}

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]
		tcur = (*pS).events[i,3]

		A = (tieval :> 0)
		if (!nowindowindeg) dynam_window_expire(pS, windeg, pwd, (*pS).gfirst[i] - 1, tcur, windowindeg, 0)
		if (!nowindowoutdeg) dynam_window_expire(pS, woutdeg, pwo, (*pS).gfirst[i] - 1, tcur, windowoutdeg, 0)
		if (nowindowindeg) {
			if (weightedindeg) indeg_i = colsum(tieval)
			else indeg_i = colsum(sign(tieval))
		}
		else {
			if (weightedindeg) indeg_i = colsum(windeg)
			else indeg_i = colsum(sign(windeg))
		}
		if (nowindowoutdeg) {
			if (weightedoutdeg) outdeg_i = rowsum(tieval)'
			else outdeg_i = rowsum(sign(tieval))'
		}
		else {
			if (weightedoutdeg) outdeg_i = rowsum(woutdeg)'
			else outdeg_i = rowsum(sign(woutdeg))'
		}
		if (active[4]) {
			AA = A * A
			nodetrans_i = rowsum(A :* AA)'
		}
		if (active[5]) {
			tertiusnumer = tertiusvec * A
			// Always the count of DISTINCT in-neighbors, never the
			// weighted indegree - see the choice engine's own identical
			// comment above dynam_choice_loglik_grad_multi().
			tertiusdenom = colsum(A)
			hasin = (tertiusdenom :> 0)
			tertius_i = hasin :* (tertiusnumer :/ (tertiusdenom :+ (1 :- hasin)))
		}

		lp = J(1, n, 0)
		if (active[1]) lp = lp :+ theta[idx[1]] :* indeg_i
		if (active[2]) lp = lp :+ theta[idx[2]] :* outdeg_i
		if (active[3]) lp = lp :+ theta[idx[3]] :* egovec
		if (active[4]) lp = lp :+ theta[idx[4]] :* nodetrans_i
		if (active[5]) lp = lp :+ theta[idx[5]] :* tertius_i
		// Two-mode: mode-2 actors (indices n1+1..n) can never be the
		// next sender - goldfish's own architecture assumes mode 1
		// only ever sends (see nwdynam.ado's own two-mode comment for
		// the real-goldfish verification this mirrors).
		if (n1 > 0) lp[(n1+1)..n] = J(1, n-n1, -1e300)

		mx = max(lp)
		P = exp(lp :- mx)
		lrsum = mx + ln(sum(P))
		ll = ll + (lp[s] - lrsum)

		P = P :/ sum(P)
		if (active[1]) gradfull[1] = gradfull[1] + (indeg_i[s] - sum(indeg_i :* P))
		if (active[2]) gradfull[2] = gradfull[2] + (outdeg_i[s] - sum(outdeg_i :* P))
		if (active[3]) gradfull[3] = gradfull[3] + (egovec[s] - sum(egovec :* P))
		if (active[4]) gradfull[4] = gradfull[4] + (nodetrans_i[s] - sum(nodetrans_i :* P))
		if (active[5]) gradfull[5] = gradfull[5] + (tertius_i[s] - sum(tertius_i :* P))

		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			a = (*pS).events[k,1]
			b = (*pS).events[k,2]
			v = (*pS).events[k,4]
			tieval[a,b] = tieval[a,b] + v
			if (!nowindowindeg) windeg[a,b] = windeg[a,b] + v
			if (!nowindowoutdeg) woutdeg[a,b] = woutdeg[a,b] + v
		}
	}

	grad = J(1, 0, 0)
	for (e=1; e<=5; e++) if (active[e]) grad = (grad, gradfull[e])
}

real scalar dynam_rate_loglik_multi(real rowvector theta, real rowvector active,
		real rowvector egovec, real rowvector tertiusvec, real scalar windowindeg, real scalar windowoutdeg,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1,
		pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_rate_loglik_grad_multi(theta, active, egovec, tertiusvec, windowindeg, windowoutdeg, weightedindeg, weightedoutdeg, n1, pS, ll, grad)
	return(ll)
}

/*
	WITH-INTERCEPT rate sub-model ("expansion batch 17", 2026-09-02,
	continuing the goldfish-comparison list) - the genuinely
	continuous-time competing-risks hazard variant this file's own
	header CORRECTION comment always said remained deferred: unlike the
	NO-intercept rate model (which reduces to the SAME ordinal partial
	likelihood as everything else in this package, using only the ORDER
	events happen in, never real elapsed time), giving an intercept makes
	the model genuinely sensitive to the ACTUAL SCALE of elapsed time
	between events (goldfish's own teaching1.Rmd vignette: "an intercept
	of -14 means a waiting time of 334 hours").

	FORMULA, derived from a real toy example against goldfish FIRST (an
	intercept is requested via a literal `1' as the first formula term,
	e.g. `dep ~ 1 + indeg(net)' - confirmed directly via goldfish's own
	`parseIntercept()'; `dep ~ 1' ALONE, with no other term, hits a real
	goldfish-internal bug unrelated to this package, `strsplit(): non-
	character argument' - not investigated further, out of scope, worked
	around by always fitting at least one real effect alongside the
	intercept in every crosscheck). Each actor i is an independent
	Poisson clock with hazard `exp(beta'X_i(t))'; the FIRST observed
	event contributes only `beta'X_{s_1}' (goldfish does not assume any
	risk period exists before the first observed event - confirmed
	directly: including a hazard-integral term for the very first event
	produced a real, wrong number until this was found and removed).
	Every SUBSEQUENT event k contributes
	  beta'X_{s_k} - (t_k - t_{k-1}) * sum_i exp(beta'X_i)
	summed over the full observed sequence - the standard competing-
	exponential-clocks log-likelihood (Cox-style), verified to match
	real goldfish EXACTLY (not within tolerance) on two independent toy
	fits (`rate ~ 1 + indeg`, `rate ~ 1 + outdeg`, real - not sequential -
	timestamps) - see dev/dynam_unit20_rateintercept_crosscheck.R/.do.
	The gradient (verified against numerical differentiation to ~1e-9,
	and confirmed ~0 at goldfish's own reported MLE) is the same
	"observed minus expected" score for every event but the first:
	  d logL_k/dbeta_j = X_{s_k,j} - dt_k * sum_i exp(beta'X_i) * X_{i,j}
	(k=1: just `X_{s_1,j}', no dt term), with the intercept itself
	simply X_{i,intercept}=1 for every actor.

	SCOPE: the intercept is ALWAYS included (theta[1], not an `active[]'
	toggle - matching goldfish's own "give `1' as a formula term"
	convention, not a boolean flag) alongside any subset of the SAME
	five rate-submodel effects (`indeg'/`outdeg'/`ego'/`nodetrans'/
	`tertius') the NO-intercept engine already has - reusing their own
	already-verified statistic-building code unchanged, only the
	AGGREGATION differs. `nodeTrans'/`tertius' remain UNWINDOWED and
	unweighted here too, matching their own scope limit in the
	NO-intercept engine.

	Window/`weighted=TRUE'/two-mode combined with the intercept
	("expansion batch 18", 2026-09-02) - each verified independently
	against real goldfish on a toy network with real, unevenly-spaced
	timestamps BEFORE this code was written (see dev/dynam_unit21_
	rateintercept_twomode_crosscheck.R, dynam_unit21b_..._window_...R,
	dynam_unit21c_..._weighted_...R): the SAME statistic swaps the
	no-intercept engine already uses (`indeg_i'/`outdeg_i' read from a
	real-time window or a cumulative weighted count instead of the
	full-history distinct-neighbor count; mode-2 actors masked to
	-1e300 in `lp' so their hazard underflows to ~0) turned out to
	compose correctly with the hazard-integral aggregation exactly as
	written, with no further correction needed - confirmed by exact
	agreement with real goldfish on all three, not merely assumed from
	the no-intercept case. Two-mode is verified here with `outdeg' only
	(goldfish itself rejects `indeg' - an "ego"-type effect - on a
	two-mode RATE sub-model; already established for the no-intercept
	engine, see dynam_rate_loglik_grad_multi()'s own header).
*/
void dynam_rateint_loglik_grad_multi(real rowvector theta, real rowvector active,
		real rowvector egovec, real rowvector tertiusvec, real rowvector modvec,
		pointer(class DynamState scalar) scalar pS, real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, j, e, tcur, tprev, dt, intercept
	real scalar windowindeg, windowoutdeg, weightedindeg, weightedoutdeg, n1
	real scalar nowindowindeg, nowindowoutdeg, g0, k, a, b, v, pwd, pwo
	real matrix tieval, windeg, woutdeg, A, AA
	real rowvector indeg_i, outdeg_i, nodetrans_i, tertius_i, tertiusnumer, tertiusdenom, hasin
	real rowvector lp, haz, gradfull, idx

	n = (*pS).n
	ll = 0
	gradfull = J(1, 6, 0)
	// tieval, windeg, woutdeg: tie values (summed event values) overall and within the windows - see
	// dynam_choice_loglik_grad_multi()
	tieval = J(n, n, 0)
	windeg = J(n, n, 0)
	woutdeg = J(n, n, 0)
	pwd = 1
	pwo = 1
	intercept = theta[1]

	windowindeg = modvec[1]
	windowoutdeg = modvec[2]
	weightedindeg = modvec[3]
	weightedoutdeg = modvec[4]
	n1 = modvec[5]
	nowindowindeg = (windowindeg >= 1e300)
	nowindowoutdeg = (windowoutdeg >= 1e300)

	idx = J(1, 5, 0)
	j = 1
	for (e=1; e<=5; e++) {
		if (active[e]) {
			j = j + 1
			idx[e] = j
		}
	}

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]
		tcur = (*pS).events[i,3]

		A = (tieval :> 0)
		if (!nowindowindeg) dynam_window_expire(pS, windeg, pwd, (*pS).gfirst[i] - 1, tcur, windowindeg, 0)
		if (!nowindowoutdeg) dynam_window_expire(pS, woutdeg, pwo, (*pS).gfirst[i] - 1, tcur, windowoutdeg, 0)
		if (nowindowindeg) {
			if (weightedindeg) indeg_i = colsum(tieval)
			else indeg_i = colsum(sign(tieval))
		}
		else {
			if (weightedindeg) indeg_i = colsum(windeg)
			else indeg_i = colsum(sign(windeg))
		}
		if (nowindowoutdeg) {
			if (weightedoutdeg) outdeg_i = rowsum(tieval)'
			else outdeg_i = rowsum(sign(tieval))'
		}
		else {
			if (weightedoutdeg) outdeg_i = rowsum(woutdeg)'
			else outdeg_i = rowsum(sign(woutdeg))'
		}
		if (active[4]) {
			AA = A * A
			nodetrans_i = rowsum(A :* AA)'
		}
		if (active[5]) {
			tertiusnumer = tertiusvec * A
			tertiusdenom = colsum(A)
			hasin = (tertiusdenom :> 0)
			tertius_i = hasin :* (tertiusnumer :/ (tertiusdenom :+ (1 :- hasin)))
		}

		lp = J(1, n, intercept)
		if (active[1]) lp = lp :+ theta[idx[1]] :* indeg_i
		if (active[2]) lp = lp :+ theta[idx[2]] :* outdeg_i
		if (active[3]) lp = lp :+ theta[idx[3]] :* egovec
		if (active[4]) lp = lp :+ theta[idx[4]] :* nodetrans_i
		if (active[5]) lp = lp :+ theta[idx[5]] :* tertius_i
		// Two-mode: mode-2 actors can never be the next sender - see
		// dynam_rate_loglik_grad_multi()'s own identical comment. Here
		// this also correctly zeroes their contribution to the hazard
		// SUM below (exp(-1e300) underflows to 0), not just the
		// softmax the no-intercept engine uses.
		if (n1 > 0) lp[(n1+1)..n] = J(1, n-n1, -1e300)

		if (i == 1) {
			ll = ll + lp[s]
			gradfull[1] = gradfull[1] + 1
			if (active[1]) gradfull[idx[1]] = gradfull[idx[1]] + indeg_i[s]
			if (active[2]) gradfull[idx[2]] = gradfull[idx[2]] + outdeg_i[s]
			if (active[3]) gradfull[idx[3]] = gradfull[idx[3]] + egovec[s]
			if (active[4]) gradfull[idx[4]] = gradfull[idx[4]] + nodetrans_i[s]
			if (active[5]) gradfull[idx[5]] = gradfull[idx[5]] + tertius_i[s]
		}
		else {
			dt = tcur - tprev
			haz = exp(lp)
			ll = ll + lp[s] - dt * sum(haz)
			gradfull[1] = gradfull[1] + (1 - dt * sum(haz))
			if (active[1]) gradfull[idx[1]] = gradfull[idx[1]] + (indeg_i[s] - dt * sum(haz :* indeg_i))
			if (active[2]) gradfull[idx[2]] = gradfull[idx[2]] + (outdeg_i[s] - dt * sum(haz :* outdeg_i))
			if (active[3]) gradfull[idx[3]] = gradfull[idx[3]] + (egovec[s] - dt * sum(haz :* egovec))
			if (active[4]) gradfull[idx[4]] = gradfull[idx[4]] + (nodetrans_i[s] - dt * sum(haz :* nodetrans_i))
			if (active[5]) gradfull[idx[5]] = gradfull[idx[5]] + (tertius_i[s] - dt * sum(haz :* tertius_i))
		}

		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			a = (*pS).events[k,1]
			b = (*pS).events[k,2]
			v = (*pS).events[k,4]
			tieval[a,b] = tieval[a,b] + v
			if (!nowindowindeg) windeg[a,b] = windeg[a,b] + v
			if (!nowindowoutdeg) woutdeg[a,b] = woutdeg[a,b] + v
		}
		tprev = tcur
	}

	grad = gradfull[1..1]
	for (e=1; e<=5; e++) if (active[e]) grad = (grad, gradfull[idx[e]])
}

real scalar dynam_rateint_loglik_multi(real rowvector theta, real rowvector active,
		real rowvector egovec, real rowvector tertiusvec, real rowvector modvec,
		pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_rateint_loglik_grad_multi(theta, active, egovec, tertiusvec, modvec, pS, ll, grad)
	return(ll)
}

/*
	CACHED evaluator of the with-intercept rate sub-model (nwdynam 1.4.0). Every statistic of the likelihood
	(indeg, outdeg, ego, nodeTrans, tertius, with windows, weights and the two-mode mask) depends only on the
	event history, never on theta. dynam_rateint_build_cache() therefore computes them ONCE per event row
	(rows of a group of equal times share the state from before the group) into N x n matrices, and
	dynam_rateint_cached_eval() evaluates log likelihood, gradient and Hessian of
	  sum_i [ lp(i, s_i) - w_i * sum_j exp(lp(i, j)) ],  lp(i, j) = theta_0 + sum_k theta_k X_k(i, j),
	w_1 = 0 (the first event has no hazard integral), w_i = t_i - t_(i-1) (so 0 within a group of equal
	times), with a few matrix operations. The function is concave in theta, so Newton-Raphson from a data-driven
	start finds the unique maximum; the analytic Hessian gives the standard errors. The formulas are those of
	dynam_rateint_loglik_grad_multi() above (the former evaluator, kept as the fallback for very large
	problems), which was verified against goldfish.
*/
struct DynamRICache {
	pointer(real matrix) rowvector X   // one N x n matrix per active effect, in the order of the active effects
	real rowvector sXs                 // per active effect: sum over events of the statistic of the sender
	real colvector w                   // time gap before each event (0 for the first event and within groups)
	real scalar N, n, n1, P
}

struct DynamRICache scalar dynam_rateint_build_cache(pointer(class DynamState scalar) scalar pS, real rowvector active,
		real rowvector egovec, real rowvector tertiusvec, real rowvector modvec) {
	struct DynamRICache scalar C
	real scalar i, s, tcur, tprev, g0, k, a, b, v, pwd, pwo, e, N, n, p
	real scalar windowindeg, windowoutdeg, weightedindeg, weightedoutdeg, nowindowindeg, nowindowoutdeg
	real matrix tieval, windeg, woutdeg, A, AA
	real rowvector indeg_i, outdeg_i, nodetrans_i, tertius_i, tertiusnumer, tertiusdenom, hasin, row
	real rowvector eff, indegrun, outdegrun

	n = (*pS).n
	N = (*pS).nevents
	C.N = N
	C.n = n
	C.n1 = modvec[5]
	eff = selectindex(active)
	p = cols(eff)
	C.P = p
	C.X = J(1, p, NULL)
	for (k=1; k<=p; k++) C.X[k] = &(J(N, n, 0))
	C.sXs = J(1, p, 0)
	C.w = J(N, 1, 0)

	tieval = J(n, n, 0)
	windeg = J(n, n, 0)
	woutdeg = J(n, n, 0)
	indegrun = J(1, n, 0)
	outdegrun = J(1, n, 0)
	pwd = 1
	pwo = 1
	windowindeg = modvec[1]
	windowoutdeg = modvec[2]
	weightedindeg = modvec[3]
	weightedoutdeg = modvec[4]
	nowindowindeg = (windowindeg >= 1e300)
	nowindowoutdeg = (windowoutdeg >= 1e300)

	for (i=1; i<=N; i++) {
		s = (*pS).events[i,1]
		tcur = (*pS).events[i,3]
		if ((*pS).gfirst[i] == i) {
			// first row of a group (or a sequential row): build the statistics of the current state
			if (active[4] | active[5]) A = (tieval :> 0)
			if (!nowindowindeg) dynam_window_expire(pS, windeg, pwd, (*pS).gfirst[i] - 1, tcur, windowindeg, 0)
			if (!nowindowoutdeg) dynam_window_expire(pS, woutdeg, pwo, (*pS).gfirst[i] - 1, tcur, windowoutdeg, 0)
			// without a window the degrees are kept up to date in O(1) per tie update (below), as the native kernel does
			if (active[1]) {
				if (nowindowindeg) indeg_i = indegrun
				else {
					if (weightedindeg) indeg_i = colsum(windeg)
					else indeg_i = colsum(sign(windeg))
				}
			}
			if (active[2]) {
				if (nowindowoutdeg) outdeg_i = outdegrun
				else {
					if (weightedoutdeg) outdeg_i = rowsum(woutdeg)'
					else outdeg_i = rowsum(sign(woutdeg))'
				}
			}
			if (active[4]) {
				AA = A * A
				nodetrans_i = rowsum(A :* AA)'
			}
			if (active[5]) {
				tertiusnumer = tertiusvec * A
				tertiusdenom = colsum(A)
				hasin = (tertiusdenom :> 0)
				tertius_i = hasin :* (tertiusnumer :/ (tertiusdenom :+ (1 :- hasin)))
			}
		}
		for (k=1; k<=p; k++) {
			e = eff[k]
			if (e == 1) row = indeg_i
			else if (e == 2) row = outdeg_i
			else if (e == 3) row = egovec
			else if (e == 4) row = nodetrans_i
			else row = tertius_i
			(*C.X[k])[i,.] = row
			C.sXs[k] = C.sXs[k] + row[s]
		}
		if (i > 1) C.w[i] = tcur - tprev
		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			a = (*pS).events[k,1]
			b = (*pS).events[k,2]
			v = (*pS).events[k,4]
			if (weightedindeg) indegrun[b] = indegrun[b] + v
			else indegrun[b] = indegrun[b] - sign(tieval[a,b]) + sign(tieval[a,b] + v)
			if (weightedoutdeg) outdegrun[a] = outdegrun[a] + v
			else outdegrun[a] = outdegrun[a] - sign(tieval[a,b]) + sign(tieval[a,b] + v)
			tieval[a,b] = tieval[a,b] + v
			if (!nowindowindeg) windeg[a,b] = windeg[a,b] + v
			if (!nowindowoutdeg) woutdeg[a,b] = woutdeg[a,b] + v
		}
		tprev = tcur
	}
	return(C)
}

void dynam_rateint_cached_eval(real scalar todo, real rowvector theta, pointer(struct DynamRICache scalar) scalar pC,
		transmorphic pAux, real scalar y, real rowvector g, real matrix H) {
	real scalar k, l, m, P, N, n, nn, sw
	real matrix L, Ew, EX
	real colvector sr
	real rowvector gg

	P = (*pC).P
	N = (*pC).N
	n = (*pC).n
	nn = ((*pC).n1 > 0 ? (*pC).n1 : n)
	L = J(N, nn, theta[1])
	for (k=1; k<=P; k++) L = L :+ theta[k+1] :* (*(*pC).X[k])[., 1..nn]
	// overflow guard: exp() of a large number is missing in Mata and sum() would drop it silently; report a missing
	// value instead so that the optimizer takes a shorter step
	m = max(L)
	if (m > 700 | missing(m)) {
		y = .
		if (todo >= 1) g = J(1, P+1, .)
		if (todo >= 2) H = J(P+1, P+1, .)
		return
	}
	// mode-2 actors of a two-mode network are not at risk: their columns are left out (nn < n)
	Ew = exp(L) :* (*pC).w
	sr = rowsum(Ew)
	sw = sum(sr)
	y = N * theta[1] + sum(theta[2..(P+1)] :* (*pC).sXs) - sw
	if (todo >= 1) {
		gg = J(1, P+1, 0)
		gg[1] = N - sw
		for (k=1; k<=P; k++) gg[k+1] = (*pC).sXs[k] - sum(Ew :* (*(*pC).X[k])[., 1..nn])
		g = gg
	}
	if (todo >= 2) {
		H = J(P+1, P+1, 0)
		H[1,1] = -sw
		for (k=1; k<=P; k++) {
			EX = Ew :* (*(*pC).X[k])[., 1..nn]
			H[1,k+1] = -sum(EX)
			H[k+1,1] = H[1,k+1]
			for (l=k; l<=P; l++) {
				H[k+1,l+1] = -sum(EX :* (*(*pC).X[l])[., 1..nn])
				H[l+1,k+1] = H[k+1,l+1]
			}
		}
	}
}

/*
	choice_coordination sub-model ("expansion batch 11", 2026-09-02, user:
	"do choice_coordination now"): goldfish's own third DyNAM sub-model, a
	genuinely different likelihood from `choice' (a single conditional
	logit over one sender's own candidate receivers) - a "multinomial-
	multinomial" joint model over UNORDERED actor pairs, for undirected
	("coordination") tie-formation events where neither actor plays a
	privileged sender/receiver role (Stadtfeld, Hollway & Block 2017,
	"Dynamic Network Actor Models: Investigating Coordination Ties through
	Time," Sociological Methodology 47(1), equations 6/9/11).

	FORMULA, derived from goldfish's own C++ source
	(src/compute_coordination_selection.cpp in the CRAN 1.6.12 source
	tarball - the paper's own prose was not relied on, the actual
	production code was read directly), then independently verified
	against a REAL goldfish fit on two toy examples (not merely read and
	trusted): for each event, build the full n x n matrix of linear
	predictors beta'S(a,b) for every ordered candidate pair a!=b (the SAME
	per-effect statistic formulas as the plain `choice' sub-model - every
	`update_DyNAM_choice_coordination_X' function in goldfish's own
	namespace is confirmed, by direct inspection, to be a THIN WRAPPER
	calling the identical `update_DyNAM_choice_X' used by plain choice;
	only the AGGREGATION differs, not the statistics). Then:
	  p(a,b) = exp(beta'S(a,b)) / sum_{b'!=a} exp(beta'S(a,b'))   -- row-
	    normalized: an ordinary choice-submodel probability that a, as
	    sender, would pick b.
	  P(a,b) = p(a,b) * p(b,a)                                    -- the
	    (unnormalized) joint mass that a and b select EACH OTHER.
	  Pnorm(a,b) = P(a,b) / (sum(P)/2)                            -- P is
	    symmetric by construction (P(a,b)=P(b,a)), so accu(P) double-
	    counts every unordered pair; dividing by 2 normalizes Pnorm to a
	    genuine probability distribution over the n*(n-1)/2 unordered
	    candidate pairs.
	  loglik contribution = ln(Pnorm(sender,receiver)).
	Verified EXACTLY (to 6 decimal places) against a real goldfish fit on
	a 4-actor/4-event toy network (inertia only, converged coefficient
	-4.172212, goldfish's own reported logLik -5.729604) and, separately,
	a 6-actor/12-event toy network with TWO effects (inertia + same(),
	converged coefficients -0.4636/-0.5638, goldfish's own reported logLik
	-31.33132) - see dev/dynam_unit14_coordination_crosscheck.R/.do.

	GRADIENT: derived by hand from the loglik formula above (not
	transcribed from goldfish's own C++ loop, though independently
	confirmed to reduce to the same quantity), then verified against
	numerical differentiation to ~3.5e-9. For each active effect k with
	per-event statistic matrix S_k(a,b):
	  dlogp_k(a,b) = S_k(a,b) - sum_{b'} p(a,b') * S_k(a,b')   -- the
	    usual "observed minus row-expected" softmax score, computed for
	    EVERY row a (not just the observed sender), since the aggregate
	    term below needs every pair's own score.
	  d loglik_event/dbeta_k = [dlogp_k(sender,receiver) +
	    dlogp_k(receiver,sender)] - sum_{a,b} Pnorm(a,b) * dlogp_k(a,b)
	The second (subtracted) term was originally derived, matching
	goldfish's own C++ structure exactly, as a sum over unordered pairs
	a<b of Pnorm(a,b)*[dlogp_k(a,b)+dlogp_k(b,a)] - algebraically
	simplified (confirmed both symbolically and via a direct numeric
	check) to the single full-matrix elementwise sum sum(Pnorm :* S_k')
	shown above, avoiding an O(n^2) Mata for-loop per effect per event.

	SCOPE (disclosed rather than silently narrow): ten effects -
	`inertia', `indeg' (alter-type, matching choice's own), `same',
	`diff', `sim', `alter' (batch 11), `nodeTrans' and `trans' (batch
	12), plus `tertius' and `four' (batch 13, added 2026-09-02) - chosen
	because each is confirmed (by direct inspection of goldfish's own
	registered effect functions, not the vignette's own
	effect-eligibility table, which self-contradicts on `alter' - the
	table lists it "x" for choice_coordination while its OWN
	accompanying prose text says the opposite, and a direct
	`preprocessingOnly=TRUE' test against real goldfish confirms the
	PROSE is correct) to be BOTH real-goldfish-eligible for
	choice_coordination AND already buildable from an existing full n x n
	state matrix (`lastcontact', or `A'/`AA'/`AtA' for the closure
	family) or trivial per-actor broadcast in this engine. `recip',
	`outdeg', `ego', `cycle', `commonSender', `commonReceiver' are
	confirmed REJECTED by real goldfish itself for choice_coordination
	(structurally sensible - an undirected tie has no direction to
	reciprocate, and indeg/outdeg coincide for an undirected network so
	only one is offered). `nodeTrans'/`tertius' turned out to need NO new
	full-matrix logic beyond what `indeg' already computes: both are
	properties of the CANDIDATE alone (`rowsum(A:*AA)'' for nodeTrans,
	the same alter-type mean-covariate ratio plain choice's own tertius
	uses), broadcasting identically. `trans' - unlike plain choice's own
	`trans_i = AA[s,.]' (one row, the actual sender) - needs the FULL
	`AA' matrix directly (`AA[a,b]' already equals the two-path count
	for ANY a as sender, b as candidate). `four' needed the genuinely
	NEW full-matrix adaptation plain choice's own distinct-actors
	correction implies: `four_i(j) = (A*(A'*A))[i,j] - A[i,j]*
	(rowsum(A[i,.]) + diagonal(A'*A)[j] - A[i,j])' generalizes to
	`S_four(a,b) = AAtA[a,b] - A[a,b]*(rowsum(A[a,.]) +
	diagonal(AtA)[b] - A[a,b])' - note the diagonal correction term is
	indexed by the CANDIDATE b, not the sender a, a genuine
	transcription pitfall caught directly (a first R draft indexed it by
	a instead, reproducing a real but WRONG number, -32.039 vs.
	goldfish's own -32.288, until corrected). All four verified directly
	against real goldfish on the SAME 6-actor/12-event toy network batch
	11 used (see dev/dynam_unit15_coordination_closure_crosscheck.R/.do
	for nodeTrans/trans, dev/dynam_unit16_coordination_tertius_four_
	crosscheck.R/.do for tertius/four). `egoAlterInt', `tie',
	`tertiusDiff', `mixedTrans' remain confirmed real-goldfish-ELIGIBLE
	(via the same direct test) but NOT YET WIRED - a real, disclosed
	follow-on for a later batch. Two-mode (bipartite) coordination
	networks are likewise out of scope - goldfish's own C++ code has a
	`twomode_or_reflexive' flag suggesting real support exists there, but
	its own joint P=p.*p' construction requires the actor-1 and actor-2
	candidate pools to be the SAME SIZE (a square matrix), a real
	structural question not yet investigated for this package's own
	two-mode convention. No windowing or `weighted=TRUE' support yet
	either (goldfish's own effect signatures allow both on several of
	these effects) - all real, disclosed gaps.
*/
void dynam_coord_loglik_grad_multi(real rowvector theta, real rowvector active,
		real matrix covmat, real matrix tiemat,
		pointer(class DynamState scalar) scalar pS, real scalar ll, real rowvector grad) {
	real scalar i, s, r, n, j, e, tcur, Z, g0, k, a, b, v
	real matrix tieval, statmat, praw, p, P, Pnorm, S_inertia, S_indeg, S_same, S_diff, S_sim, S_alter
	real matrix A, AA, AtA, AAtA, S_nodetrans, S_trans, S_tertius, S_four, S_egoalterint, S_tie
	real rowvector samevec, diffvec, simvec, altervec, tertiusvec, egoalterint1vec, egoalterint2vec
	real rowvector indegvec, nodetransvec, idx, gradfull
	real rowvector tertiusnumer, tertiusdenom, hasin
	real colvector rmax, rs, dg
	real matrix dlogp, colrep, rowrep

	samevec = covmat[1,.]
	diffvec = covmat[2,.]
	simvec = covmat[3,.]
	altervec = covmat[4,.]
	tertiusvec = covmat[5,.]
	egoalterint1vec = covmat[6,.]
	egoalterint2vec = covmat[7,.]

	n = (*pS).n
	ll = 0
	gradfull = J(1, 12, 0)
	// tieval: tie values of the undirected network (summed event values, both cells); a tie exists while above 0
	tieval = J(n, n, 0)

	idx = J(1, 12, 0)
	j = 0
	for (e=1; e<=12; e++) {
		if (active[e]) {
			j = j + 1
			idx[e] = j
		}
	}

	// Mata's `:' elementwise operators do NOT broadcast a plain row
	// vector against its own transpose (1 x n vs n x 1) the way they
	// broadcast an n x n matrix against an n x 1 column (confirmed
	// directly - a first draft assumed the former worked too and hit a
	// real conformability error) - so the outer same/diff/sim/alter
	// matrices are built explicitly via J()-replication instead.
	// same()/diff()/sim()/alter() never depend on event history, so
	// they are computed ONCE here, outside the event loop, not
	// recomputed every event.
	colrep = samevec' * J(1, n, 1)
	rowrep = J(n, 1, 1) * samevec
	S_same = (colrep :== rowrep)
	colrep = diffvec' * J(1, n, 1)
	rowrep = J(n, 1, 1) * diffvec
	S_diff = abs(colrep :- rowrep)
	colrep = simvec' * J(1, n, 1)
	rowrep = J(n, 1, 1) * simvec
	S_sim = -abs(colrep :- rowrep)
	S_alter = J(n, 1, 1) * altervec
	// egoAlterInt: sender a's own covariate1 TIMES candidate b's own
	// covariate2 - an outer product, matching plain choice's own
	// `egoalterint_i = egoalterint1vec[s] :* egoalterint2vec' exactly
	// (that IS the outer product's row s already). Static like
	// same/diff/sim/alter, computed once here.
	S_egoalterint = egoalterint1vec' * egoalterint2vec
	// tie() ("cross-network effects, v1 scope", extended to
	// choice_coordination 2026-09-02): presence in a SEPARATE, STATIC
	// exogenous network - identical in shape to plain choice's own
	// `tie_i = tiemat[s,.]' but, since choice_coordination already needs
	// the FULL n x n statistic matrix for every effect (not just one
	// row), `S_tie' is simply `tiemat' itself - no computation at all,
	// the simplest possible addition. `tiemat' has 0 rows when tie() was
	// not given (matching every other "not given" sentinel convention in
	// this engine). Verified directly against real goldfish (see
	// dev/dynam_unit19_coordination_tie_crosscheck.R/.do).
	if (active[12]) S_tie = tiemat

	for (i=1; i<=(*pS).nevents; i++) {
		s = (*pS).events[i,1]
		r = (*pS).events[i,2]
		tcur = (*pS).events[i,3]

		S_inertia = (tieval :> 0)
		indegvec = colsum(sign(tieval))
		S_indeg = J(n, 1, 1) * indegvec
		if (active[7] | active[8] | active[9] | active[10]) {
			A = S_inertia
			AA = A * A
		}
		// nodeTrans (candidate b's own embeddedness in transitive
		// structures) - a property of b alone, not a sender/candidate
		// comparison, exactly like `indeg' - broadcasts across every row
		// a. `trans' (two-paths sender->k->candidate) - unlike the plain
		// choice engine's own `trans_i = AA[s,.]' (one row, the actual
		// sender), the FULL matrix version needed here is simply AA
		// itself: AA[a,b] already equals the two-path count for ANY a
		// as sender, b as candidate. Both verified directly against real
		// goldfish (dev/dynam_unit15_coordination_closure_crosscheck.R/
		// .do) before being wired in.
		if (active[7]) {
			nodetransvec = rowsum(A :* AA)'
			S_nodetrans = J(n, 1, 1) * nodetransvec
		}
		if (active[8]) S_trans = AA
		// tertius (candidate b's own mean covariate over ITS in-
		// neighbors) - like nodeTrans, a property of b alone (matching
		// plain choice's own tertius_i, which likewise never varies by
		// sender), so broadcasts across every row a exactly like
		// nodeTrans/indeg. Imputed 0 when b has no in-neighbors yet,
		// matching plain choice's own verified-not-documented behavior.
		if (active[9]) {
			tertiusnumer = tertiusvec * A
			tertiusdenom = colsum(A)
			hasin = (tertiusdenom :> 0)
			S_tertius = J(n, 1, 1) * (hasin :* (tertiusnumer :/ (tertiusdenom :+ (1 :- hasin))))
		}
		// four (three-paths a->k<-l->b, k!=l, all four actors distinct)
		// - unlike nodeTrans/tertius, this genuinely needs the FULL
		// matrix, since (per plain choice's own already-verified
		// distinct-actors correction) the subtracted term depends on
		// BOTH a (`rowsum(A[a,.])') and b (`diagonal(AtA)[b]') - the
		// diagonal term is indexed by the CANDIDATE b, not the sender a
		// (confirmed directly: an initial R draft indexed it by a
		// instead by mistranscribing plain choice's own `diagonal(AtA)
		// [j]' where j is the candidate, and it reproduced a real but
		// WRONG number until corrected - see dev/dynam_unit16_
		// coordination_tertius_four_crosscheck.R/.do).
		if (active[10]) {
			AtA = A' * A
			AAtA = A * AtA
			rs = rowsum(A)
			dg = diagonal(AtA)
			S_four = AAtA :- A :* ((rs * J(1,n,1)) :+ (J(n,1,1) * dg') :- A)
		}

		statmat = J(n, n, 0)
		if (active[1]) statmat = statmat :+ theta[idx[1]] :* S_inertia
		if (active[2]) statmat = statmat :+ theta[idx[2]] :* S_indeg
		if (active[3]) statmat = statmat :+ theta[idx[3]] :* S_same
		if (active[4]) statmat = statmat :+ theta[idx[4]] :* S_diff
		if (active[5]) statmat = statmat :+ theta[idx[5]] :* S_sim
		if (active[6]) statmat = statmat :+ theta[idx[6]] :* S_alter
		if (active[7]) statmat = statmat :+ theta[idx[7]] :* S_nodetrans
		if (active[8]) statmat = statmat :+ theta[idx[8]] :* S_trans
		if (active[9]) statmat = statmat :+ theta[idx[9]] :* S_tertius
		if (active[10]) statmat = statmat :+ theta[idx[10]] :* S_four
		if (active[11]) statmat = statmat :+ theta[idx[11]] :* S_egoalterint
		if (active[12]) statmat = statmat :+ theta[idx[12]] :* S_tie

		rmax = rowmax(statmat)
		praw = exp(statmat :- rmax)
		_diag(praw, 0)
		p = praw :/ rowsum(praw)

		P = p :* p'
		Z = sum(P) / 2
		Pnorm = P :/ Z

		ll = ll + ln(Pnorm[s,r])

		if (active[1]) {
			dlogp = S_inertia :- rowsum(p :* S_inertia)
			gradfull[1] = gradfull[1] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[2]) {
			dlogp = S_indeg :- rowsum(p :* S_indeg)
			gradfull[2] = gradfull[2] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[3]) {
			dlogp = S_same :- rowsum(p :* S_same)
			gradfull[3] = gradfull[3] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[4]) {
			dlogp = S_diff :- rowsum(p :* S_diff)
			gradfull[4] = gradfull[4] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[5]) {
			dlogp = S_sim :- rowsum(p :* S_sim)
			gradfull[5] = gradfull[5] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[6]) {
			dlogp = S_alter :- rowsum(p :* S_alter)
			gradfull[6] = gradfull[6] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[7]) {
			dlogp = S_nodetrans :- rowsum(p :* S_nodetrans)
			gradfull[7] = gradfull[7] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[8]) {
			dlogp = S_trans :- rowsum(p :* S_trans)
			gradfull[8] = gradfull[8] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[9]) {
			dlogp = S_tertius :- rowsum(p :* S_tertius)
			gradfull[9] = gradfull[9] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[10]) {
			dlogp = S_four :- rowsum(p :* S_four)
			gradfull[10] = gradfull[10] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[11]) {
			dlogp = S_egoalterint :- rowsum(p :* S_egoalterint)
			gradfull[11] = gradfull[11] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}
		if (active[12]) {
			dlogp = S_tie :- rowsum(p :* S_tie)
			gradfull[12] = gradfull[12] + ((dlogp[s,r] + dlogp[r,s]) - sum(Pnorm :* dlogp))
		}

		if ((g0 = dynam_group_end(pS, i)) > 0) for (k=g0; k<=i; k++) {
			a = (*pS).events[k,1]
			b = (*pS).events[k,2]
			v = (*pS).events[k,4]
			tieval[a,b] = tieval[a,b] + v
			tieval[b,a] = tieval[b,a] + v
		}
	}

	grad = J(1, 0, 0)
	for (e=1; e<=12; e++) if (active[e]) grad = (grad, gradfull[e])
}

real scalar dynam_coord_loglik_multi(real rowvector theta, real rowvector active,
		real matrix covmat, real matrix tiemat,
		pointer(class DynamState scalar) scalar pS) {
	real scalar ll
	real rowvector grad
	dynam_coord_loglik_grad_multi(theta, active, covmat, tiemat, pS, ll, grad)
	return(ll)
}

void dynam_choice_eval_unit1(real scalar todo, real rowvector theta, pointer(class DynamState scalar) scalar pS,
		real scalar y, real rowvector g, real matrix H) {
	real rowvector grad
	dynam_choice_loglik_grad_unit1(theta, pS, y, grad)
	if (todo >= 1) g = grad
}

void dynam_rate_eval_unit1(real scalar todo, real rowvector theta, pointer(class DynamState scalar) scalar pS,
		real scalar y, real rowvector g, real matrix H) {
	real rowvector grad
	dynam_rate_loglik_grad_unit1(theta, pS, y, grad)
	if (todo >= 1) g = grad
}

end
