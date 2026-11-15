/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
 * nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	unw_dynam.do -- Mata part of nwdynam (Dynamic Network Actor Models, derived from the R package goldfish).

	their effects, windows, weighted effects, two-mode networks, opportunities() and tie(). This file prepares the event stream
	(DynamState: events sorted by time, simultaneous events grouped as in goldfish), writes the model into the plugin's data
	frame (class DynamNat), runs Mata's optimize() on the plugin evaluator and returns the estimates to the ado file. There is no
	Mata fallback: when the plugin is missing or too old the command stops (r(601), r(3495)).

	The Mata evaluators that carried the likelihoods before (verified against goldfish; the effect definitions and their history are
	documented there) are reference code of the test suite only: cscripts/ref_dynam_mata.do. The plugin is tested against them.
*/

capture mata: mata drop DynamState()
capture mata: mata drop dynam_node_rows()
capture mata: mata drop DynamNativePluginFilename()
capture mata: mata drop DynamNativePluginSubdir()
capture mata: mata drop DynamNativePluginPath()
capture mata: mata drop DynamNativeAvailable()
capture mata: mata drop DynamNat()
capture mata: mata drop DYNAM_PROTOCOL_NEEDED()
capture mata: mata drop dynam_nat_eval()
capture mata: mata drop dynam_nat_eval_nr()
capture mata: mata drop dynam_nat_eval_d1()
capture mata: mata drop dynam_nat_optimize()
capture mata: mata drop dynam_rateint_accept()
capture mata: mata drop dynam_rateint_newton()
capture mata: mata drop DynamFitUnit1()
capture mata: mata drop DynamFitRateUnit1()
capture mata: mata drop DynamChoiceFitMulti()
capture mata: mata drop DynamRateFitMulti()
capture mata: mata drop DynamRateInterceptFitMulti()




capture mata: mata drop dynam_group_end()
capture mata: mata drop dynam_coord_loglik_grad_multi()
capture mata: mata drop dynam_coord_loglik_multi()
capture mata: mata drop dynam_coord_eval_multi()
capture mata: mata drop DynamCoordFitMulti()

mata:
mata set matastrict off


/*
	DynamState -- the sorted event stream. Deliberately simpler than
	unw_rem.do's own RemState (no cumulative degree accumulator
	matrices precomputed up front) - the choice sub-model's own three
	v1 effects are all functions of a single incrementally-maintained
	BINARY tie matrix (see dynam_choice_loglik_grad_unit1() below), not
	of raw per-actor event counts the way REM's NODSnd/NIDRec are, so
	there is nothing analogous to build_degree_accumulators() to
	precompute here.
*/
class DynamState {
	real matrix events     // nevents x 4: (sender, receiver, time, value), SORTED ascending by time; value 1 when the
	                       // network has no event values (as goldfish, the value is the increment of the tie)
	real scalar n          // number of actors
	real scalar nevents
	real colvector gfirst  // for each event: row of the first event of its group of equal times (its own row if sequential)
	real scalar sequential // 1 = events with equal times are processed one after another (old behaviour)
	real scalar ngroups    // number of distinct time groups
	real scalar nativevalued // 1 = the native kernel gets the event values (v5): some value is 0 or below

	void init()
}

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
	dynam_node_rows(): for every node name, its row in the label column lab of the data (0 when absent, -1 when the
	label occurs more than once). Used by nwdynam.ado to match per-actor covariates to the actors by label.
*/
real colvector dynam_node_rows(string rowvector names, string colvector lab)
{
	transmorphic A
	real colvector pos
	real scalar k

	A = asarray_create("string", 1)
	for (k = 1; k <= rows(lab); k++) {
		if (asarray_contains(A, lab[k])) asarray(A, lab[k], -1)
		else asarray(A, lab[k], k)
	}
	pos = J(cols(names), 1, 0)
	for (k = 1; k <= cols(names); k++) {
		if (asarray_contains(A, names[k])) pos[k] = asarray(A, names[k])
	}
	return(pos)
}

void DynamState::init(real matrix rawevents, real scalar nn) {
	real colvector ord
	real scalar k

	// Matching unw_rem.do's own RemState::init() self-loop rejection -
	// DyNAM's choice sub-model risk set is defined over i != j
	// candidate receivers only (a sender cannot choose itself), so a
	// self-loop event has no well-defined receiver-set position and
	// must be rejected outright, not silently dropped.
	for (k=1; k<=rows(rawevents); k++) {
		if (rawevents[k,1] == rawevents[k,2]) {
			// the message goes through errprintf(): _error() fails with r(3300) on texts longer than about 110 characters
			errprintf("{p 0 0 2}nwdynam: event %s has sender == receiver (self-loop); DyNAM's choice sub-model is defined over i != j dyads only.{p_end}\n", strofreal(k))
			exit(3498)
		}
	}

	n = nn
	// get_eventlist() makes no ordering guarantee (see unw_rem.do's own
	// RemState::init() comment, verified there by reading
	// unw_core.do's nwattime_slice_event()) - the choice sub-model's
	// own incrementally-updated tie matrix requires chronological
	// order.
	// Events with the SAME time keep the order in which the network lists them: order() on the time column alone breaks
	// ties arbitrarily (it depends on the state of the random-number generator), so the same data gave different tie
	// histories and a different likelihood from run to run. The row number as second key makes the sort stable.
	ord = order((rawevents[.,3], (1::rows(rawevents))), (1, 2))
	events = rawevents[ord,.]
	nevents = rows(events)
	// column 4: the value of each event (an event list without values, from an older core, counts every event as 1)
	if (cols(events) < 4) events = (events[., 1..3], J(nevents, 1, 1))
	else events = events[., 1..4]
	// option sequential (scalar __nwdynam_sequential = 1): old behaviour, every event is its own group
	sequential = (st_numscalar("__nwdynam_sequential") == 1)
	gfirst = (1::nevents)
	if (!sequential) {
		for (k=2; k<=nevents; k++) {
			if (events[k,3] == events[k-1,3]) gfirst[k] = gfirst[k-1]
		}
	}
	ngroups = sum(gfirst :== (1::nevents))
	nativevalued = 0
}

/*
	Native (C) plugin dispatch (native/dynam_sim.c) - own dedicated
	functions, mirroring unw_ergm.do's own ErgmNative* pattern (a single
	plugin serving one command's own engine) rather than unw_core.do's
	shared NativeGraph* dispatcher (which serves several graph-algorithm
	commands sharing one plugin) - DyNAM is its own third/fourth engine,
	same "why a separate engine" reasoning the development notes already
	gives for keeping unw_dynam.do independent of unw_rem.do.

	UNLIKE ergm_mcmc.c/saom_sim.c (a full sampling/simulation loop lives
	entirely in C, called once per outer MCMC-step/ministep from Mata),
	dynam_sim.c accelerates only the log-likelihood/gradient EVALUATION
	at a given theta - Mata's own optimize() still drives the outer
	BFGS/Nelder-Mead loop unchanged (see dynam_sim.c's own header for
	why: reimplementing an optimizer in C is a real, avoidable risk this
	design does not need to take). `DynamNativeAvailable()` below
	follows the exact same "return 0, never error, caller falls back to
	Mata" contract as `ErgmNativeAvailable()`/`NativeGraphAvailable()`.

	TWO genuinely different lookup strategies are needed (harmonisation
	2026-09-02, the development notes), tried in order:
	 (1) `findfile()` on the platform-specific basename alone - what
	     actually finds the plugin after a real `net install`, which
	     flattens every package "f" line into
	     PLUS/<firstletter-of-basename>/<basename>, discarding any
	     declared subdirectory entirely (verified directly by installing
	     a real test package into a scratch PLUS directory). Distinct
	     per-platform basenames (macOS and Windows used to share the
	     bare "dynam_sim.plugin" name; only Unix had its own "_unix"
	     suffix) are exactly what let this same flat PLUS folder hold
	     all three platforms' binaries at once without collision.
	 (2) a manually-constructed path relative to nwdynam.ado's own
	     directory, `lib/plugins/<os>/<name>' - unreachable after a real
	     net install (per (1) above, since net install has no
	     subdirectory concept at all) but still needed for a raw git
	     checkout (`adopath ++ <repo-root>`, this project's own
	     dev-mode/regression-testing convention): `findfile()` does not
	     search subdirectories of a plain adopath entry (confirmed
	     directly - only PLUS's own single-letter-subfolder convention
	     gets that treatment), so the nested lib/plugins/<os>/ layout
	     the repo itself uses is invisible to strategy (1) alone.
*/
string scalar DynamNativePluginFilename(){
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("dynam_sim_windows.plugin")
	if (os == "Unix") return("dynam_sim_unix.plugin")
	return("dynam_sim_macos.plugin")
}

string scalar DynamNativePluginSubdir(){
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("windows")
	if (os == "Unix") return("unix")
	return("macos")
}

string scalar DynamNativePluginPath(){
	string scalar fname, found, full, dir, fn

	fname = DynamNativePluginFilename()
	found = findfile(fname)
	if (found != "") return(found)

	full = findfile("nwdynam.ado")
	if (full == "") return("")
	pathsplit(full, dir, fn)
	return(pathjoin(pathjoin(dir, "lib"),
		pathjoin("plugins", pathjoin(DynamNativePluginSubdir(), fname))))
}

real scalar DynamNativeAvailable(){
	string scalar p

	p = DynamNativePluginPath()
	if (p == "") return(0)
	return(fileexists(p))
}

/*
	Acceptance test for a Newton-Raphson fit of the with-intercept rate sub-model: at the optimum the Newton decrement
	g' (-H)^-1 g is zero, and optimize() may stop only after tolerances that round-off cannot meet at large log
	likelihoods (it then keeps taking "backed up" steps). A result is accepted when the decrement is below 1e-8, i.e.
	the log likelihood is within 1e-8 of its maximum.
*/
real scalar dynam_rateint_accept(transmorphic S_opt, real scalar errcode) {
	real rowvector g
	real matrix H, Hi

	if (errcode != 0 & errcode != 430) return(0)
	g = optimize_result_gradient(S_opt)
	H = optimize_result_Hessian(S_opt)
	if (missing(g) | missing(H)) return(0)
	Hi = invsym(-H)
	if (min(diagonal(Hi)) <= 0) return(0)
	return(g * Hi * g' < 1e-8)
}

/*
	Newton-Raphson on a concave log likelihood with analytic gradient and Hessian: the evaluator f has two argument
	slots a1, a2. optimize() gets close (and stops after a Newton decrement below 1e-11), but it only accepts a step
	that raises the log likelihood, which round-off defeats a few digits before the maximum (about 1e-8 in the
	parameters at a log likelihood of -100, worse at -400000). Explicit Newton steps, judged by the size of the gradient
	instead, then bring the gradient to round-off level. ok = 1 when a result was found; theta_hat, V (the inverse of
	minus the Hessian, i.e. the observed information as optimize_result_V_oim() would give it) and ll are then set.
*/
void dynam_rateint_newton(transmorphic f, real rowvector theta0, transmorphic a1, transmorphic a2,
		real rowvector theta_hat, real matrix V, real scalar ll, real scalar ok) {
	transmorphic S_opt
	real rowvector t, tn, g, gn
	real matrix H, Hn
	real scalar errcode, y, yn, it, converged

	ok = 0
	S_opt = optimize_init()
	optimize_init_evaluator(S_opt, f)
	optimize_init_evaluatortype(S_opt, "d2")
	optimize_init_params(S_opt, theta0)
	optimize_init_argument(S_opt, 1, a1)
	optimize_init_argument(S_opt, 2, a2)
	optimize_init_technique(S_opt, "nr")
	optimize_init_which(S_opt, "max")
	optimize_init_tracelevel(S_opt, "none")
	optimize_init_conv_warning(S_opt, "off")
	optimize_init_conv_maxiter(S_opt, 40)
	optimize_init_conv_ptol(S_opt, 1e-13)
	optimize_init_conv_vtol(S_opt, 1e-14)
	optimize_init_conv_nrtol(S_opt, 1e-11)
	errcode = _optimize(S_opt)
	if (!dynam_rateint_accept(S_opt, errcode)) return
	t = optimize_result_params(S_opt)
	(*f)(2, t, a1, a2, y, g, H)
	if (missing(y) | missing(g) | missing(H)) return
	for (it=1; it<=5; it++) {
		tn = t + g * invsym(-H)
		(*f)(2, tn, a1, a2, yn, gn, Hn)
		if (missing(yn) | missing(gn) | missing(Hn)) break
		if (sqrt(gn * gn') >= sqrt(g * g')) break
		t = tn
		y = yn
		g = gn
		H = Hn
	}
	V = invsym(-H)
	if (min(diagonal(V)) <= 0) return
	theta_hat = t
	ll = y
	ok = 1
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

void dynam_coord_eval_multi(real scalar todo, real rowvector theta, pointer(real rowvector) scalar pActive,
		pointer(real matrix) scalar pCovmat, pointer(real matrix) scalar pTiemat,
		pointer(class DynamState scalar) scalar pS, real scalar y, real rowvector g, real matrix H) {
	real rowvector grad
	dynam_coord_loglik_grad_multi(theta, *pActive, *pCovmat, *pTiemat, pS, y, grad)
	if (todo >= 1) g = grad
}

/*
	DynamCoordFitMulti() -- fits any nonempty subset of choice_coordination's
	eight effects (six from "batch 11," plus `nodeTrans'/`trans' added
	"batch 12" - see dynam_coord_loglik_grad_multi()'s own header
	comment for the formula, its real-goldfish verification, and the
	disclosed remaining scope). `eventmat' rows are unordered actor pairs (any
	order - the likelihood is symmetric in sender/receiver by
	construction, confirmed directly), one per undirected tie-formation
	event, matching `nwset ..., eventtime() undirected''s own raw
	`get_eventlist()' output (confirmed directly: no auto-symmetrization,
	one row per declared event). Same random-perturbation multi-start
	strategy as `DynamChoiceFitMulti'/`DynamRateFitMulti'. Only 3
	`optimize_init_argument()' slots used (activevec, covmat, S) - no
	risk yet of the 9-slot cap that forced the `covmat'/`windowvec'/
	`weightvec' bundling in the choice engine, but the SAME bundling
	convention (one matrix argument for per-actor covariates) is used
	from the start here too, so future growth (windowing, weighting,
	closure effects) never needs a signature change.
*/
void DynamCoordFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector samevec, real rowvector diffvec, real rowvector simvec, real rowvector altervec,
		real rowvector tertiusvec, real rowvector egoalterint1vec, real rowvector egoalterint2vec,
		real matrix tiemat,
		string scalar bname, string scalar vname, string scalar llname) {
	class DynamState scalar S
	transmorphic S_opt
	real rowvector theta0, theta_hat
	real matrix V_best, covmat
	real scalar nparams, attempt, nconverged, errcode, nstarts, thisval, bestval
	real matrix starts

	S = DynamState()
	S.init(eventmat, n)

	covmat = (samevec \ diffvec \ simvec \ altervec \ tertiusvec \ egoalterint1vec \ egoalterint2vec)

	nparams = sum(activevec)
	if (nparams == 0) _error("nwdynam: at least one effect must be selected.")

	nstarts = 8
	starts = J(nstarts, nparams, 0)
	for (attempt=1; attempt<=nstarts; attempt++) {
		starts[attempt,.] = (runiform(1,nparams) :- 0.5) :* (0.15 * (mod(attempt-1, nstarts) + 1))
	}

	// Unlike the plain choice/rate sub-models (an ordinary softmax,
	// globally concave, so the FIRST attempt _optimize() completes
	// without error is trustworthy as the unique optimum), the
	// coordination sub-model's own P=p.*p' joint construction is a
	// PRODUCT of two softmax probabilities, not guaranteed concave -
	// found directly, not assumed: a first draft here accepted any
	// `_optimize()' call that returned errcode 0, which does NOT mean
	// the optimizer actually reached a stationary point (Stata's own
	// `_optimize()' returns 0 whenever the numerical procedure
	// completes without a fatal error, including stopping early at
	// `conv_maxiter' with a real gradient still far from zero) - this
	// silently accepted a converged-looking but WRONG two-effect
	// coordinate on a real toy example (inertia+same: landed at
	// (-0.399,0.056), true joint MLE (-0.464,-0.564), confirmed by
	// hand-evaluating this file's own dynam_coord_loglik_grad_multi()
	// directly at both points - the reported gradient at the accepted
	// point was (0.04,-3.64), nowhere near zero). Fixed two ways: (1)
	// require `optimize_result_converged()' == 1, not just errcode==0;
	// (2) run EVERY start (not just until the first "success") and
	// keep the BEST (highest loglik) among the genuinely converged
	// ones, since a non-concave surface can have more than one
	// stationary point and the first one reached is not guaranteed
	// global. A SECOND real bug found while fixing the first: keeping
	// only a saved `S_opt' HANDLE from the winning attempt and reading
	// `optimize_result_V_oim()'/`optimize_result_value()' off it AFTER
	// the loop ended silently returned a DIFFERENT (later) attempt's
	// own results instead - confirmed directly by tracing the starts
	// matrix and each attempt's own converged value, which showed the
	// function's own final answer matching a LATER attempt's results,
	// not the one flagged as best. `optimize_init()' handles are not
	// safe to keep across later `optimize_init()' calls reusing the
	// same Mata variable - every value needed from the winning attempt
	// (`theta_hat', `V_best', `bestval') is now extracted into a plain
	// value INSIDE the loop, the instant that attempt is found to be
	// the new best, never deferred to after the loop.
	nconverged = 0
	bestval = .
	for (attempt=1; attempt<=rows(starts); attempt++) {
		theta0 = starts[attempt,.]
		S_opt = optimize_init()
		optimize_init_evaluator(S_opt, &dynam_coord_eval_multi())
		optimize_init_evaluatortype(S_opt, "d1")
		optimize_init_params(S_opt, theta0)
		optimize_init_argument(S_opt, 1, &activevec)
		optimize_init_argument(S_opt, 2, &covmat)
		optimize_init_argument(S_opt, 3, &tiemat)
		optimize_init_argument(S_opt, 4, &S)
		optimize_init_technique(S_opt, "bfgs")
		optimize_init_which(S_opt, "max")
		optimize_init_tracelevel(S_opt, "none")
		optimize_init_conv_warning(S_opt, "off")
		optimize_init_conv_maxiter(S_opt, 100)
		errcode = _optimize(S_opt)
		if (errcode == 0 & optimize_result_converged(S_opt) == 1) {
			thisval = optimize_result_value(S_opt)
			if (nconverged == 0 | thisval > bestval) {
				bestval = thisval
				theta_hat = optimize_result_params(S_opt)
				V_best = optimize_result_V_oim(S_opt)
			}
			nconverged = nconverged + 1
		}
	}
	if (nconverged == 0) {
		S_opt = optimize_init()
		optimize_init_evaluator(S_opt, &dynam_coord_eval_multi())
		optimize_init_evaluatortype(S_opt, "d1")
		optimize_init_params(S_opt, J(1, nparams, 0.01))
		optimize_init_argument(S_opt, 1, &activevec)
		optimize_init_argument(S_opt, 2, &covmat)
		optimize_init_argument(S_opt, 3, &tiemat)
		optimize_init_argument(S_opt, 4, &S)
		optimize_init_technique(S_opt, "nm")
		optimize_init_which(S_opt, "max")
		optimize_init_tracelevel(S_opt, "none")
		optimize_init_conv_warning(S_opt, "off")
		optimize_init_conv_maxiter(S_opt, 300)
		theta_hat = optimize(S_opt)
		bestval = optimize_result_value(S_opt)
		V_best = optimize_result_V_oim(S_opt)
	}

	st_matrix(bname, theta_hat)
	st_matrix(vname, V_best)
	st_numscalar(llname, bestval)
}
/*
	NATIVE ENGINE. Every nwdynam fit runs in the compiled plugin (native/dynam_sim.c); there is no Mata fallback.
	class DynamNat writes the model (events, covariates, windows, masks) once per fit into the native frame
	__nwdynam_native (variable v1, "blob"; layout in the header of dynam_sim.c) and dynam_nat_eval() is the optimize()
	evaluator: one plugin call per evaluation, only theta travels (at full double precision), the result (log
	likelihood, gradient, Hessian) comes back in v2. Mata's optimize() drives the outer loop.

	The plugin must exist and speak the protocol of this package (DynamNat::open() asks it with code 99):
	  r(601)  the plugin file for this platform is not installed
	  r(3495) the installed plugin is too old for this version of nwdynam (or cannot be loaded)
	Both stop the command with a message; the Mata evaluators of the earlier versions survive only as reference code of
	the test suite (cscripts/ref_dynam_mata.do).
*/
class DynamNat {
	real scalar L        // length of the blob
	real scalar Pt       // number of parameters
	real scalar isopen
	string scalar origframe

	void open()
	void close()
	void setup()
}

real scalar DYNAM_PROTOCOL_NEEDED() return(2)

void DynamNat::open() {
	string scalar p, fname, q
	real scalar rc, junk, ver

	isopen = 0
	fname = DynamNativePluginFilename()
	p = DynamNativePluginPath()
	if (p == "" | !fileexists(p)) {
		errprintf("{p 0 0 2}nwdynam: the compiled plugin {bf:%s} was not found. nwdynam fits every model with this plugin and has no Mata fallback; reinstall the package with {bf:net install nwdynam, replace}.{p_end}\n", fname)
		exit(601)
	}
	origframe = st_framecurrent()
	stata("capture frame drop __nwdynam_native")
	stata("frame create __nwdynam_native")
	st_framecurrent("__nwdynam_native")
	isopen = 1
	st_addobs(4)
	junk = st_addvar("double", "v1")
	junk = st_addvar("double", "v2")
	q = char(34)
	stata("capture program drop dynamnative")
	stata("capture program dynamnative, plugin using(" + q + p + q + ")")
	// the zeros after the code are for plugins of version 1, which read three numbers before they look at the code (and crash
	// on a shorter string); they answer 198 to code 99
	rc = _stata("plugin call dynamnative v1 v2, " + q + "99 0 0 0 0 0 0" + q, 1)
	ver = .
	if (rc == 0) ver = st_data(1, 2)
	if (rc != 0 | ver != DYNAM_PROTOCOL_NEEDED()) {
		close()
		errprintf("{p 0 0 2}nwdynam: the installed plugin {bf:%s} is too old for this version of nwdynam (or it cannot be loaded). Reinstall the package with {bf:net install nwdynam, replace}.{p_end}\n", fname)
		exit(3495)
	}
}

void DynamNat::close() {
	if (isopen) {
		st_framecurrent(origframe)
		stata("capture frame drop __nwdynam_native")
		isopen = 0
	}
}

/*
	setup(): model 1 choice, 2 rate, 3 rate with intercept, 4 coordination. active = 0/1 vector of the effect slots,
	covmat = 7 x n covariate rows, win4/wgt4 = windows and weighted flags (inertia, recip, indeg, outdeg), n1 = mode-1
	actors (0 = one-mode), oppmat (nevents x n or 0 x 0), tiemat (n x n or 0 x 0).
*/
void DynamNat::setup(real scalar model, pointer(class DynamState scalar) scalar pS, real rowvector active,
		real matrix covmat, real rowvector win4, real rowvector wgt4, real scalar n1, real matrix oppmat, real matrix tiemat) {
	real scalar nev, hasopp, hastie, M, nres
	real colvector gend, blob
	real rowvector hdr

	nev = (*pS).nevents
	M = cols(active)
	hasopp = (rows(oppmat) > 0)
	hastie = (rows(tiemat) > 0)
	gend = J(nev, 1, 1)
	if (nev > 1) gend[1::(nev-1)] = ((*pS).gfirst[2::nev] :== (2::nev))
	hdr = (DYNAM_PROTOCOL_NEEDED(), model, (*pS).n, nev, n1, M, hasopp, hastie, win4, wgt4, active, J(1, 20 - M, 0))
	blob = hdr' \ vec(((*pS).events, gend)') \ vec(covmat')
	if (hasopp) blob = blob \ vec(oppmat')
	if (hastie) blob = blob \ vec(tiemat')
	L = rows(blob)
	Pt = sum(active) + (model == 3)
	nres = 1 + Pt + Pt * Pt
	if (max((L, nres)) > st_nobs()) st_addobs(max((L, nres)) - st_nobs())
	st_store((1::L), 1, blob)
}

void dynam_nat_eval(real scalar todo, real rowvector theta, pointer(class DynamNat scalar) scalar pN,
		real scalar y, real rowvector g, real matrix H) {
	real scalar k, pt
	string scalar argstr
	real colvector res

	pt = (*pN).Pt
	argstr = "20 " + strofreal((*pN).L, "%18.0f") + " " + strofreal(min((todo, 2))) + " " + strofreal(pt)
	for (k = 1; k <= pt; k++) argstr = argstr + " " + strofreal(theta[k], "%25.17g")
	stata("plugin call dynamnative v1 v2, " + char(34) + argstr + char(34))
	res = st_data((1::(1 + pt + (todo >= 2 ? pt * pt : 0))), 2)
	y = res[1]
	if (missing(y)) {
		if (todo >= 1) g = J(1, pt, .)
		if (todo >= 2) H = J(pt, pt, .)
		return
	}
	if (todo >= 1) g = res[2::(1 + pt)]'
	if (todo >= 2) H = rowshape(res[(2 + pt)::(1 + pt + pt * pt)], pt)
}

// the same evaluator with the two argument slots of dynam_rateint_newton()
void dynam_nat_eval_nr(real scalar todo, real rowvector theta, pointer(class DynamNat scalar) scalar pN, transmorphic pAux,
		real scalar y, real rowvector g, real matrix H) {
	dynam_nat_eval(todo, theta, pN, y, g, H)
}

void dynam_nat_eval_d1(real scalar todo, real rowvector theta, pointer(class DynamNat scalar) scalar pN,
		real scalar y, real rowvector g, real matrix H) {
	dynam_nat_eval(min((todo, 1)), theta, pN, y, g, H)
}

/*
	Maximisation on the native evaluator. mode 0 (choice and rate, concave): up to 8 random starts, BFGS with the analytic
	Hessian, the first converged start wins, Nelder-Mead from 0.01 as the last resort. mode 1 (coordination, not concave):
	every start is run and the best converged result is kept (see the reference evaluator for the history of that rule).
	mode 2 (rate with intercept, concave): Newton-Raphson from a data-driven start, then BFGS from it, then mode 0.
	Fills theta_hat, V, ll.
*/
void dynam_nat_optimize(pointer(class DynamNat scalar) scalar pN, real scalar mode, real rowvector theta0data,
		real rowvector theta_hat, real matrix V, real scalar ll) {
	transmorphic S_opt
	real rowvector theta0, th2
	real matrix starts, V_best, V2
	real scalar nparams, nstarts, attempt, ok, errcode, nconverged, thisval, bestval, nrok, ll2

	nparams = (*pN).Pt
	nstarts = 8
	if (mode == 2) {
		dynam_rateint_newton(&dynam_nat_eval_nr(), theta0data, pN, 0, theta_hat, V, ll, nrok)
		if (nrok) return
		S_opt = optimize_init()
		optimize_init_evaluator(S_opt, &dynam_nat_eval())
		optimize_init_evaluatortype(S_opt, "d2")
		optimize_init_params(S_opt, theta0data)
		optimize_init_argument(S_opt, 1, pN)
		optimize_init_technique(S_opt, "bfgs")
		optimize_init_which(S_opt, "max")
		optimize_init_tracelevel(S_opt, "none")
		optimize_init_conv_warning(S_opt, "off")
		optimize_init_conv_maxiter(S_opt, 200)
		errcode = _optimize(S_opt)
		if (errcode == 0 & optimize_result_converged(S_opt) == 1) {
			theta_hat = optimize_result_params(S_opt)
			if (!missing(theta_hat)) {
				V = optimize_result_V_oim(S_opt)
				if (!missing(V)) {
					ll = optimize_result_value(S_opt)
					return
				}
			}
		}
		mode = 0
	}
	starts = J(nstarts, nparams, 0)
	for (attempt = 1; attempt <= nstarts; attempt++) {
		starts[attempt, .] = (runiform(1, nparams) :- 0.5) :* (0.15 * (mod(attempt - 1, nstarts) + 1))
	}
	// the first start is fixed, so that a fit without seed() is reproducible
	starts[1, .] = J(1, nparams, 0.01)
	nconverged = 0
	bestval = .
	ok = 0
	for (attempt = 1; attempt <= nstarts & (mode == 1 | ok == 0); attempt++) {
		theta0 = starts[attempt, .]
		S_opt = optimize_init()
		if (mode == 1) optimize_init_evaluator(S_opt, &dynam_nat_eval_d1())
		else optimize_init_evaluator(S_opt, &dynam_nat_eval())
		optimize_init_evaluatortype(S_opt, (mode == 1 ? "d1" : "d2"))
		optimize_init_params(S_opt, theta0)
		optimize_init_argument(S_opt, 1, pN)
		optimize_init_technique(S_opt, "bfgs")
		optimize_init_which(S_opt, "max")
		optimize_init_tracelevel(S_opt, "none")
		optimize_init_conv_warning(S_opt, "off")
		optimize_init_conv_maxiter(S_opt, (mode == 1 ? 100 : 60))
		errcode = _optimize(S_opt)
		if (mode == 1) {
			if (errcode == 0 & optimize_result_converged(S_opt) == 1) {
				thisval = optimize_result_value(S_opt)
				if (nconverged == 0 | thisval > bestval) {
					bestval = thisval
					theta_hat = optimize_result_params(S_opt)
					V_best = optimize_result_V_oim(S_opt)
				}
				nconverged = nconverged + 1
			}
		}
		else if (errcode == 0) {
			theta_hat = optimize_result_params(S_opt)
			V = optimize_result_V_oim(S_opt)
			ll = optimize_result_value(S_opt)
			ok = 1
		}
	}
	if (mode == 1 & nconverged > 0) {
		V = V_best
		ll = bestval
		return
	}
	if (ok) {
		// concave log likelihood: Newton-Raphson on the analytic Hessian from the BFGS result brings the gradient to round-off
		// level (BFGS stops at its default tolerances); the BFGS result stays when Newton does not accept
		dynam_rateint_newton(&dynam_nat_eval_nr(), theta_hat, pN, 0, th2, V2, ll2, nrok)
		if (nrok) {
			theta_hat = th2
			V = V2
			ll = ll2
		}
		return
	}
	S_opt = optimize_init()
	optimize_init_evaluator(S_opt, &dynam_nat_eval_d1())
	optimize_init_evaluatortype(S_opt, "d1")
	optimize_init_params(S_opt, J(1, nparams, 0.01))
	optimize_init_argument(S_opt, 1, pN)
	optimize_init_technique(S_opt, "nm")
	optimize_init_which(S_opt, "max")
	optimize_init_tracelevel(S_opt, "none")
	optimize_init_conv_warning(S_opt, "off")
	optimize_init_conv_maxiter(S_opt, 300)
	theta_hat = optimize(S_opt)
	V = optimize_result_V_oim(S_opt)
	ll = optimize_result_value(S_opt)
}

/* ---- the fit functions called by nwdynam.ado (signatures unchanged since 1.4.0) ---- */

void DynamChoiceFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector samevec, real rowvector diffvec, real rowvector simvec, real rowvector altervec,
		real rowvector tertiusvec, real rowvector egoalterint1vec, real rowvector egoalterint2vec,
		real scalar windowinertia, real scalar windowrecip,
		real scalar weightedinertia, real scalar weightedrecip, real scalar weightedindeg, real scalar weightedoutdeg,
		real scalar n1, real matrix oppmat, real matrix tiemat,
		string scalar bname, string scalar vname, string scalar llname) {
	class DynamState scalar S
	class DynamNat scalar NT
	real matrix covmat, V
	real rowvector theta_hat
	real scalar ll

	S = DynamState()
	S.init(eventmat, n)
	if (sum(activevec) == 0) _error("nwdynam: at least one effect must be selected.")
	covmat = (samevec \ diffvec \ simvec \ altervec \ tertiusvec \ egoalterint1vec \ egoalterint2vec)
	NT = DynamNat()
	NT.open()
	NT.setup(1, &S, activevec, covmat, (windowinertia, windowrecip, 1e300, 1e300),
		(weightedinertia, weightedrecip, weightedindeg, weightedoutdeg), n1, oppmat, tiemat)
	dynam_nat_optimize(&NT, 0, ., theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

void DynamRateFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector egovec, real rowvector tertiusvec, real scalar windowindeg, real scalar windowoutdeg,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1,
		string scalar bname, string scalar vname, string scalar llname) {
	class DynamState scalar S
	class DynamNat scalar NT
	real matrix covmat, V
	real rowvector theta_hat
	real scalar ll

	S = DynamState()
	S.init(eventmat, n)
	if (sum(activevec) == 0) _error("nwdynam: at least one effect must be selected.")
	covmat = (J(3, n, 0) \ egovec \ tertiusvec \ J(2, n, 0))
	NT = DynamNat()
	NT.open()
	NT.setup(2, &S, activevec, covmat, (1e300, 1e300, windowindeg, windowoutdeg),
		(0, 0, weightedindeg, weightedoutdeg), n1, J(0, 0, 0), J(0, 0, 0))
	dynam_nat_optimize(&NT, 0, ., theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

void DynamRateInterceptFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector egovec, real rowvector tertiusvec,
		real scalar windowindeg, real scalar windowoutdeg,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1,
		string scalar bname, string scalar vname, string scalar llname) {
	class DynamState scalar S
	class DynamNat scalar NT
	real matrix covmat, V
	real rowvector theta_hat, theta0
	real scalar ll, nev, ttot, nrisk

	S = DynamState()
	S.init(eventmat, n)
	covmat = (J(3, n, 0) \ egovec \ tertiusvec \ J(2, n, 0))
	nev = S.nevents
	ttot = S.events[nev, 3] - S.events[1, 3]
	nrisk = (n1 > 0 ? n1 : n)
	// start: the intercept of the model without effects (events per actor and time unit), all effects 0
	theta0 = J(1, 1 + sum(activevec), 0)
	if (ttot > 0 & nev > 1) theta0[1] = ln((nev - 1) / (nrisk * ttot))
	NT = DynamNat()
	NT.open()
	NT.setup(3, &S, activevec, covmat, (1e300, 1e300, windowindeg, windowoutdeg),
		(0, 0, weightedindeg, weightedoutdeg), n1, J(0, 0, 0), J(0, 0, 0))
	dynam_nat_optimize(&NT, 2, theta0, theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

// the three-effect choice model and the two-effect rate model without options (kept as short forms)
void DynamFitUnit1(real matrix eventmat, real scalar n, string scalar bname, string scalar vname, string scalar llname) {
	DynamChoiceFitMulti(eventmat, n, (1, 1, 1, J(1, 14, 0)), J(1, n, 0), J(1, n, 0), J(1, n, 0), J(1, n, 0), J(1, n, 0), J(1, n, 0), J(1, n, 0),
		1e300, 1e300, 0, 0, 0, 0, 0, J(0, 0, 0), J(0, 0, 0), bname, vname, llname)
}

void DynamFitRateUnit1(real matrix eventmat, real scalar n, string scalar bname, string scalar vname, string scalar llname) {
	DynamRateFitMulti(eventmat, n, (1, 1, 0, 0, 0), J(1, n, 0), J(1, n, 0), 1e300, 1e300, 0, 0, 0, bname, vname, llname)
}

end
