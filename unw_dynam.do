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

	Every fit runs in the compiled plugin native/dynam_sim.c: choice, rate, rate with intercept and choice_coordination, with
	their effects, windows, weighted effects, two-mode networks, opportunities() and tie(). This file prepares the event stream
	(DynamState: events sorted by time, simultaneous events grouped as in goldfish), writes the model into the plugin's data
	frame (class DynamNat), runs Mata's optimize() on the plugin evaluator and returns the estimates to the ado file. There is no
	Mata fallback: when the plugin is missing or too old the command stops (r(601), r(3495)).

	The Mata evaluators that carried the likelihoods before (verified against goldfish; the effect definitions and their history are
	documented there) are reference code of the test suite only: cscripts/ref_dynam_mata.do. The plugin is tested against them.
*/

capture mata: mata drop DynamState()
capture mata: mata drop dynam_node_rows()
capture mata: mata drop dynam_change_rows()
capture mata: mata drop dynam_comp_rows()
capture mata: mata drop dynam_check_presence()
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

/*
	dynam_change_rows(): the change events of one time-varying covariate as rows (covariate row, actor, time, new value), sorted by
	time (the row order of the event list breaks ties). ev = event list of the change network (id, id, time, value) with the
	actor ids of that network, chgnames = its node names, mainnames = the node names of the fitted network, init = the values
	at the start (one per actor of the fitted network), slots = the covariate rows (1..7) that use the variable. With
	increment = 1 the values are increments of the current value (goldfish: increment), else new values (replace). Stops with
	r(198) for an event whose sender is not its receiver, an actor that is not an actor of the fitted network or a missing value.
*/
real matrix dynam_change_rows(real matrix ev, string rowvector chgnames, string rowvector mainnames, real rowvector init,
		real scalar increment, real rowvector slots, string scalar netname, string scalar varname) {
	real colvector ord, pos
	real matrix e, out
	real rowvector cur
	real scalar k, a, v, q
	transmorphic A

	if (rows(ev) == 0) return(J(0, 4, 0))
	if (cols(ev) < 4) {
		errprintf("{p 0 0 2}nwdynam: change(): %s has no event values; declare it with a value variable ({bf:nwset id id value, eventtime(time)}): the value is the new value of %s (or its increment with option increment).{p_end}\n", netname, varname)
		exit(198)
	}
	A = asarray_create("string", 1)
	for (k = 1; k <= cols(mainnames); k++) asarray(A, mainnames[k], k)
	e = ev[., 1..4]
	ord = order((e[., 3], (1::rows(e))), (1, 2))
	e = e[ord, .]
	pos = J(rows(e), 1, 0)
	for (k = 1; k <= rows(e); k++) {
		if (e[k, 1] != e[k, 2]) {
			errprintf("{p 0 0 2}nwdynam: change(): an event of %s has sender %s and receiver %s; a change event names one actor (sender = receiver), for example {bf:nwset id id value, eventtime(time) selfloop}.{p_end}\n", netname, chgnames[e[k, 1]], chgnames[e[k, 2]])
			exit(198)
		}
		if (missing(e[k, 4])) {
			errprintf("{p 0 0 2}nwdynam: change(): an event of %s has a missing value.{p_end}\n", netname)
			exit(198)
		}
		if (!asarray_contains(A, chgnames[e[k, 1]])) {
			errprintf("{p 0 0 2}nwdynam: change(): actor %s of %s is not an actor of the fitted network (actors are matched by label).{p_end}\n", chgnames[e[k, 1]], netname)
			exit(198)
		}
		pos[k] = asarray(A, chgnames[e[k, 1]])
	}
	cur = init
	out = J(rows(e) * cols(slots), 4, .)
	for (k = 1; k <= rows(e); k++) {
		a = pos[k]
		v = (increment ? cur[a] + e[k, 4] : e[k, 4])
		cur[a] = v
		for (q = 1; q <= cols(slots); q++) out[(k - 1) * cols(slots) + q, .] = (slots[q], a, e[k, 3], v)
	}
	return(out)
}

/*
	dynam_comp_rows(): composition change. ev = event list of the composition network (id, id, time, value) with the actor ids of that
	network, chgnames = its node names, mainnames = node names of the fitted network, init = presence at the start (one 0/1 value
	per actor of the fitted network; empty = everybody is present), netname for messages. The value of an event is the new presence
	of the actor (1 = joins, 0 = leaves; goldfish: replace = TRUE / FALSE); sender = receiver. Returns rows (8, actor, time, value)
	sorted by time (the row order of the event list breaks ties) and, when init is given, rows (9, actor, ., value) for every actor.
	Stops with r(198) for an event whose sender is not its receiver, an actor that is not an actor of the fitted network, a missing
	value or a value other than 0 and 1.
*/
real matrix dynam_comp_rows(real matrix ev, string rowvector chgnames, string rowvector mainnames, real rowvector init,
		string scalar netname) {
	real colvector ord, pos
	real matrix e, out
	real scalar k, n
	transmorphic A

	n = cols(mainnames)
	out = J(0, 4, 0)
	if (cols(init) > 0) out = (J(n, 1, 9), (1::n), J(n, 1, .), init')
	if (rows(ev) == 0) return(out)
	if (cols(ev) < 4) {
		errprintf("{p 0 0 2}nwdynam: compchange(): %s has no event values; declare it with a value variable ({bf:nwset id id value, eventtime(time) selfloop}): 1 = the actor joins, 0 = the actor leaves.{p_end}\n", netname)
		exit(198)
	}
	A = asarray_create("string", 1)
	for (k = 1; k <= n; k++) asarray(A, mainnames[k], k)
	e = ev[., 1..4]
	ord = order((e[., 3], (1::rows(e))), (1, 2))
	e = e[ord, .]
	pos = J(rows(e), 1, 0)
	for (k = 1; k <= rows(e); k++) {
		if (e[k, 1] != e[k, 2]) {
			errprintf("{p 0 0 2}nwdynam: compchange(): an event of %s has sender %s and receiver %s; a composition change names one actor (sender = receiver), for example {bf:nwset id id value, eventtime(time) selfloop}.{p_end}\n", netname, chgnames[e[k, 1]], chgnames[e[k, 2]])
			exit(198)
		}
		if (missing(e[k, 4]) | (e[k, 4] != 0 & e[k, 4] != 1)) {
			errprintf("{p 0 0 2}nwdynam: compchange(): an event of %s has the value %s; the value must be 1 (the actor joins) or 0 (the actor leaves).{p_end}\n", netname, strofreal(e[k, 4]))
			exit(198)
		}
		if (!asarray_contains(A, chgnames[e[k, 1]])) {
			errprintf("{p 0 0 2}nwdynam: compchange(): actor %s of %s is not an actor of the fitted network (actors are matched by label).{p_end}\n", chgnames[e[k, 1]], netname)
			exit(198)
		}
		pos[k] = asarray(A, chgnames[e[k, 1]])
	}
	return(out \ (J(rows(e), 1, 8), pos, e[., 3], e[., 4]))
}

/*
	dynam_check_presence(): as goldfish, an event whose sender or receiver is absent at the time of the event is an error (the
	change events at the time of the event count: presence = the last change at or before that time). rows = the 8/9 rows of
	dynam_comp_rows(), events = the event list (sender, receiver, time, ...), names = actor names of the fitted network, chg = the
	change events of covariates (covariate row, actor, time, value; 0 x 4 when none): goldfish refuses a change event of an absent actor too.
*/
void dynam_check_presence(real matrix rows, real matrix events, string rowvector names, real matrix chg) {
	real colvector ord, pres
	real matrix pc, pin, e
	real scalar k, p, n, q
	
	n = cols(names)
	pres = J(n, 1, 1)
	pin = select(rows, rows[., 1] :== 9)
	if (rows(pin) > 0) pres[pin[., 2]] = pin[., 4]
	pc = select(rows, rows[., 1] :== 8)
	if (rows(pc) == 0 & all(pres)) return
	e = events[., 1..3]
	ord = order((e[., 3], (1::rows(e))), (1, 2))
	e = e[ord, .]
	p = 1
	for (k = 1; k <= rows(e); k++) {
		while (p <= rows(pc)) {
			if (pc[p, 3] > e[k, 3]) break
			pres[pc[p, 2]] = pc[p, 4]
			p++
		}
		if (!pres[e[k, 1]] | !pres[e[k, 2]]) {
			errprintf("{p 0 0 2}nwdynam: actor %s is not present at time %s but has an event (sender %s, receiver %s); events of an absent actor are not allowed - add its join before the event or remove the event.{p_end}\n", names[pres[e[k, 1]] ? e[k, 2] : e[k, 1]], strofreal(e[k, 3]), names[e[k, 1]], names[e[k, 2]])
			exit(198)
		}
	}
	// the change events of covariates (sorted by time already)
	pres = J(n, 1, 1)
	if (rows(pin) > 0) pres[pin[., 2]] = pin[., 4]
	p = 1
	for (k = 1; k <= rows(chg); k++) {
		if (chg[k, 1] > 7 & chg[k, 1] != 15) continue
		while (p <= rows(pc)) {
			if (pc[p, 3] > chg[k, 3]) break
			pres[pc[p, 2]] = pc[p, 4]
			p++
		}
		if (!pres[chg[k, 2]]) {
			errprintf("{p 0 0 2}nwdynam: actor %s is not present at time %s but has a covariate change event (change()); a change event of an absent actor is not allowed - add its join before the change or remove the change.{p_end}\n", names[chg[k, 2]], strofreal(chg[k, 3]))
			exit(198)
		}
	}
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
	real rowvector t, tn, g, gn, d
	real matrix H, Hn, Hi
	real scalar y, yn, it, k, step, dec, good

	ok = 0
	t = theta0
	(*f)(2, t, a1, a2, y, g, H)
	if (missing(y) | missing(g) | missing(H)) return
	// damped Newton steps on the exact Hessian: a full step, halved until the log likelihood does not fall. The loop
	// ends when the Newton decrement g' (-H)^-1 g is below 1e-12, or 1e-13 times the size of the log likelihood where round-off
	// does not allow less (the log likelihood is then at its maximum to that precision);
	// optimize()'s own "nr" technique keeps backing up at round-off level without ever meeting its tolerances, which cost
	// up to a thousand evaluations on dense data.
	good = 0
	for (it = 1; it <= 100; it++) {
		Hi = invsym(-H)
		if (min(diagonal(Hi)) <= 0) return
		d = g * Hi
		dec = d * g'
		if (dec < max((1e-12, 1e-13 * abs(y)))) {
			good = 1
			break
		}
		step = 1
		for (k = 1; k <= 30; k++) {
			tn = t + step * d
			(*f)(2, tn, a1, a2, yn, gn, Hn)
			if (!missing(yn) & !missing(gn) & !missing(Hn)) {
				if (yn >= y) break
			}
			step = step / 2
		}
		if (k > 30) {
			// no step raises the log likelihood: only round-off is left when the decrement is already tiny
			good = (dec < 1e-8)
			break
		}
		t = tn
		y = yn
		g = gn
		H = Hn
	}
	if (!good) return
	// bring the gradient to round-off level: steps are taken while they reduce its norm
	for (it = 1; it <= 5; it++) {
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
	real matrix chg      // change events: (covariate row 1..7, or 15 = row 5 used by a windowed tertius only, actor, time, new value) of time-varying covariates, sorted by time;
	                     // rows with code 8 = composition change (8, actor, time, 1 joins / 0 leaves), code 9 = presence at the start
	                     // (9, actor, ., 0/1); unset = none

	void open()
	void close()
	void setup()
}

real scalar DYNAM_PROTOCOL_NEEDED() return(5)

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
	covmat = 7 x n covariate rows, win5 = windows (inertia, recip, indeg, outdeg, tertius), wgt4 = weighted flags (inertia, recip, indeg, outdeg), n1 = mode-1
	actors (0 = one-mode), oppmat (nevents x n or 0 x 0), tiemat (n x n or 0 x 0).
*/
void DynamNat::setup(real scalar model, pointer(class DynamState scalar) scalar pS, real rowvector active,
		real matrix covmat, real rowvector win5, real rowvector wgt4, real scalar n1, real matrix oppmat, real matrix tiemat) {
	real scalar nev, hasopp, hastie, M, nres, nchg, ncomp, hasinit
	real colvector gend, blob, pv
	real rowvector hdr
	real matrix cvc, pc, pin

	if (cols(win5) == 4) win5 = (win5, 1e300)   // the window of tertius is optional
	nev = (*pS).nevents
	M = cols(active)
	hasopp = (rows(oppmat) > 0)
	hastie = (rows(tiemat) > 0)
	gend = J(nev, 1, 1)
	if (nev > 1) gend[1::(nev-1)] = ((*pS).gfirst[2::nev] :== (2::nev))
	cvc = J(0, 4, 0); pc = J(0, 3, 0); pin = J(0, 2, 0)
	if (cols(chg) == 4 & rows(chg) > 0) {
		cvc = select(chg, (chg[., 1] :<= 7) :| (chg[., 1] :== 15))
		pc = select(chg, chg[., 1] :== 8)[., (2, 3, 4)]
		pin = select(chg, chg[., 1] :== 9)[., (2, 4)]
	}
	nchg = rows(cvc)
	ncomp = rows(pc)
	hasinit = (rows(pin) > 0)
	hdr = (DYNAM_PROTOCOL_NEEDED(), model, (*pS).n, nev, n1, M, hasopp, hastie, win5[1..4], wgt4, active, J(1, 20 - M, 0), nchg, ncomp, hasinit, win5[5])
	blob = hdr' \ vec(((*pS).events, gend)')
	if (nchg > 0) blob = blob \ vec(cvc')
	if (ncomp > 0) blob = blob \ vec(pc')
	if (hasinit) {
		// presence at the start, one value per actor (actors without a row are present)
		pv = J((*pS).n, 1, 1)
		pv[pin[., 1]] = pin[., 2]
		blob = blob \ pv
	}
	blob = blob \ vec(covmat')
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
	Maximisation on the native evaluator. mode 0 (choice and rate, concave): Newton-Raphson from 0, else up to 8 starts of BFGS
	with the analytic Hessian (the first converged start wins, then Newton polishes it), Nelder-Mead from 0.01 as the last resort. mode 1 (coordination, not concave):
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
	if (mode == 0) {
		// concave: Newton-Raphson on the exact Hessian from 0; when it is not accepted (no maximum, effects that are not
		// identified) the BFGS multi-start below runs as before
		dynam_rateint_newton(&dynam_nat_eval_nr(), J(1, nparams, 0), pN, 0, theta_hat, V, ll, nrok)
		if (nrok) return
	}
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
		real scalar windowinertia, real scalar windowrecip, real scalar windowtertius,
		real scalar weightedinertia, real scalar weightedrecip, real scalar weightedindeg, real scalar weightedoutdeg,
		real scalar n1, real matrix oppmat, real matrix tiemat, real matrix chgmat,
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
	NT.chg = chgmat
	NT.open()
	NT.setup(1, &S, activevec, covmat, (windowinertia, windowrecip, 1e300, 1e300, windowtertius),
		(weightedinertia, weightedrecip, weightedindeg, weightedoutdeg), n1, oppmat, tiemat)
	dynam_nat_optimize(&NT, 0, ., theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

void DynamRateFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector egovec, real rowvector tertiusvec, real scalar windowindeg, real scalar windowoutdeg, real scalar windowtertius,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1, real matrix chgmat,
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
	NT.chg = chgmat
	NT.open()
	NT.setup(2, &S, activevec, covmat, (1e300, 1e300, windowindeg, windowoutdeg, windowtertius),
		(0, 0, weightedindeg, weightedoutdeg), n1, J(0, 0, 0), J(0, 0, 0))
	dynam_nat_optimize(&NT, 0, ., theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

void DynamRateInterceptFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector egovec, real rowvector tertiusvec,
		real scalar windowindeg, real scalar windowoutdeg, real scalar windowtertius,
		real scalar weightedindeg, real scalar weightedoutdeg, real scalar n1, real matrix chgmat,
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
	NT.chg = chgmat
	NT.open()
	NT.setup(3, &S, activevec, covmat, (1e300, 1e300, windowindeg, windowoutdeg, windowtertius),
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
		1e300, 1e300, 1e300, 0, 0, 0, 0, 0, J(0, 0, 0), J(0, 0, 0), J(0, 4, 0), bname, vname, llname)
}

void DynamFitRateUnit1(real matrix eventmat, real scalar n, string scalar bname, string scalar vname, string scalar llname) {
	DynamRateFitMulti(eventmat, n, (1, 1, 0, 0, 0), J(1, n, 0), J(1, n, 0), 1e300, 1e300, 1e300, 0, 0, 0, J(0, 4, 0), bname, vname, llname)
}

void DynamCoordFitMulti(real matrix eventmat, real scalar n, real rowvector activevec,
		real rowvector samevec, real rowvector diffvec, real rowvector simvec, real rowvector altervec,
		real rowvector tertiusvec, real rowvector egoalterint1vec, real rowvector egoalterint2vec,
		real scalar windowinertia, real scalar windowindeg, real scalar windowtertius, real scalar weightedinertia, real scalar weightedindeg,
		real matrix tiemat, real matrix chgmat,
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
	NT.chg = chgmat
	NT.open()
	NT.setup(4, &S, activevec, covmat, (windowinertia, 1e300, windowindeg, 1e300, windowtertius), (weightedinertia, 0, weightedindeg, 0), 0, J(0, 0, 0), tiemat)
	// not concave (product of two softmax probabilities): every start is run, the best converged result is kept
	dynam_nat_optimize(&NT, 1, ., theta_hat, V, ll)
	NT.close()
	st_matrix(bname, theta_hat)
	st_matrix(vname, V)
	st_numscalar(llname, ll)
}

end
