*! version 1.9.1  18nov2026 nwdynam
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
* nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

capture program drop nwdynam
program nwdynam, eclass
	version 14
	if word(`"`0'"', 1) == "version" {
		nwdynam_version
		exit
	}
	// replay: nwdynam with no arguments, or with display options only (ci, level(), noci, nopvalues, cformat() and the
	// other options of Stata's ereturn display), after an nwdynam fit; also what estimates replay runs. Without an earlier nwdynam fit, nwdynam alone fits the current network as before.
	if replay() & "`e(cmd)'" == "nwdynam" {
		capture _nwdynam_tabopts __t check `0'
		if !_rc {
			_nwdynam_display `0'
			exit
		}
	}
	syntax [anything(name=netname)] [, SUBMODEL(string) INERTIA RECIP INDEG OUTDEG TRANS CYCLE COMMONSENDER COMMONRECEIVER FOUR NODETRANS SAME(varname numeric) DIFF(varname numeric) SIM(varname numeric) EGO(varname numeric) ALTER(varname numeric) TERTIUS(varname numeric) EGOALTERINT(varlist min=2 max=2) INERTIAWINDOW(string) RECIPWINDOW(string) INDEGWINDOW(string) OUTDEGWINDOW(string) TERTIUSWINDOW(string) WEIGHTEDINERTIA WEIGHTEDRECIP WEIGHTEDINDEG WEIGHTEDOUTDEG OPPORTUNITIES(string) TIE(string) CHange(string) INCREMENT PResent(varname numeric) COMPchange(string) OPPFRAME(string) NETTERMS(string asis) IGNOREREP(string) INTERCEPT SEED(string) SEQUENTIAL *]
	// display options of the coefficient table (ci, level(), noci, nopvalues, cformat(), ...); any other option is
	// refused here with Stata's usual message
	local __tabopts `"`options'"'
	_nwdynam_tabopts __chk check, `__tabopts'
	// window options: a positive number when given (internally -1 = not given)
	foreach __nwdynam_w in inertiawindow recipwindow indegwindow outdegwindow tertiuswindow {
		if `"``__nwdynam_w''"' != "" {
			capture confirm number ``__nwdynam_w''
			if _rc | real(`"``__nwdynam_w''"') <= 0 {
				di as error "{p 0 0 2}`__nwdynam_w'() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
				error 198
			}
		}
		else local `__nwdynam_w' = -1
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
	}
	set more off
	_nwdynam_libcheck

	_nwsyntax `netname', max(1)

	if "`egoalterint'" != "" {
		foreach __nwdynam_v of local egoalterint {
			capture confirm numeric variable `__nwdynam_v'
			if _rc {
				di as error "{p 0 0 2}egoalterint(): {bf:`__nwdynam_v'} is a string variable; egoalterint() needs two numeric variables (make a numeric variable first, for example with {help encode}).{p_end}"
				error 109
			}
		}
	}
	// present() / compchange(): composition change (actors who join or leave the network; goldfish: the node attribute present).
	// present() = the presence at the start (a 0/1 variable with one row per actor, like the covariates), compchange() = an event
	// network of the joins (value 1) and leaves (value 0), one actor per event (sender = receiver)
	if "`opportunities'" != "" & "`present'" != "" {
		di as error "{p 0 0 2}present() cannot be combined with opportunities(): both read the data in memory in a different layout - see {help nwdynam}.{p_end}"
		error 198
	}
	// opportunities() needs one row per (event, available actor); the covariate options need one row per actor
	if "`oppframe'" != "" & "`opportunities'" == "" {
		di as error "{p 0 0 2}oppframe() names the frame that holds the variables of opportunities(); give opportunities() too.{p_end}"
		error 198
	}
	if "`oppframe'" != "" {
		capture confirm frame `oppframe'
		if _rc {
			di as error "{p 0 0 2}oppframe(): frame {bf:`oppframe'} does not exist.{p_end}"
			error 111
		}
	}
	if "`opportunities'" != "" & "`oppframe'" == "" & "`same'`diff'`sim'`alter'`ego'`tertius'`egoalterint'" != "" {
		di as error "{p 0 0 2}opportunities() reads one row per event and available actor, the covariate options read one row per actor: with a covariate option, keep the covariates in memory and put the opportunities data in another frame, {bf:oppframe(}{it:framename}{bf:)} - see {help nwdynam}.{p_end}"
		error 198
	}
	if "`opportunities'" != "" {
		if `: word count `opportunities'' != 2 {
			di as error "{p 0 0 2}opportunities() takes two variables: the event sequence number and the available actor.{p_end}"
			if `: word count `opportunities'' < 2 error 102
			error 103
		}
		if "`oppframe'" == "" {
			unab opportunities : `opportunities'
			capture confirm numeric variable `opportunities'
			if _rc error 109
		}
		else {
			capture frame `oppframe': confirm numeric variable `opportunities'
		}
		if _rc {
			di as error "{p 0 0 2}opportunities(): both variables must be numeric variables of the " cond("`oppframe'" == "", "current dataset", "frame `oppframe'") ".{p_end}"
			error 111
		}
	}


	// change(): time-varying covariates, pairs "variable network": the variable is a covariate of same()/diff()/sim()/ego()/alter()/
	// tertius()/egoalterint(), the network an event network (nwset id id value, eventtime(time)) with its changes
	local __nwdynam_nchgpairs = 0
	if `"`change'"' != "" {
		local __nwdynam_nw : word count `change'
		if mod(`__nwdynam_nw', 2) != 0 {
			di as error "{p 0 0 2}change() takes pairs: a covariate variable and the event network with its changes, for example {bf:change(x xchange)}.{p_end}"
			error 198
		}
		if "`same'`diff'`sim'`ego'`alter'`tertius'`egoalterint'" == "" {
			di as error "{p 0 0 2}change() needs a covariate option (same(), diff(), sim(), ego(), alter(), tertius() or egoalterint()) that uses the variable.{p_end}"
			error 198
		}
		if "`opportunities'" != "" & "`oppframe'" == "" {
			di as error "{p 0 0 2}change() with opportunities() needs oppframe(): the covariates stay in memory, the opportunities data goes to another frame - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_nchgpairs = `__nwdynam_nw' / 2
		local __nwdynam_cvused "`same' `diff' `sim' `ego' `alter' `tertius' `egoalterint'"
		forvalues __nwdynam_i = 1/`__nwdynam_nchgpairs' {
			local __nwdynam_cvar : word `=2 * `__nwdynam_i' - 1' of `change'
			local __nwdynam_cnet : word `=2 * `__nwdynam_i'' of `change'
			capture confirm numeric variable `__nwdynam_cvar'
			if _rc {
				di as error "{p 0 0 2}change(): {bf:`__nwdynam_cvar'} is not a numeric variable of the current dataset.{p_end}"
				error 111
			}
			unab __nwdynam_cvar : `__nwdynam_cvar'
			local __nwdynam_ok : list posof "`__nwdynam_cvar'" in __nwdynam_cvused
			if `__nwdynam_ok' == 0 {
				di as error "{p 0 0 2}change(): {bf:`__nwdynam_cvar'} is not used by same(), diff(), sim(), ego(), alter(), tertius() or egoalterint(); name a variable of an effect.{p_end}"
				error 198
			}
			if `tertiuswindow' != -1 & "`intercept'" != "" & "`__nwdynam_cvar'" == "`tertius'" & "`submodel'" == "rate" {
				di as error "{p 0 0 2}change() of the tertius() variable together with tertiuswindow() and intercept in a rate model is not available: goldfish itself gives an estimate there that depends on how it orders the changes and the expiries of the window.{p_end}"
				error 198
			}
			local __nwdynam_cvars "`__nwdynam_cvars' `__nwdynam_cvar'"
			local __nwdynam_cnets "`__nwdynam_cnets' `__nwdynam_cnet'"
			local __nwdynam_dup : list dups __nwdynam_cvars
			if "`__nwdynam_dup'" != "" {
				di as error "{p 0 0 2}change(): variable {bf:`__nwdynam_dup'} is given twice; one network holds all its changes.{p_end}"
				error 198
			}
		}
	}
	else if "`increment'" != "" {
		di as error "increment needs change() - see {help nwdynam}."
		error 198
	}

	// `__nwdynam_n1' (mode-1 actor count, contiguous indices 1..n1 in
	// this package's own combined actor-index convention) - 0 sentinel
	// for "not two-mode", matching this file's own -1/1e300 "not given"
	// sentinel conventions elsewhere. Extracted unconditionally right
	// after _nwsyntax (cheap even when unused) since both the two-mode
	// effect-eligibility checks below and the Mata dispatch at the very
	// end need it.
	if "`is2mode'" == "true" {
		mata: st_numscalar("__nwdynam_n1", `netobj'->get_nodes_mode1())
		local __nwdynam_n1 = __nwdynam_n1
		scalar drop __nwdynam_n1
	}
	else {
		local __nwdynam_n1 = 0
	}

	if "`istemporal'" != "true" | "`temporaltype'" != "event" {
		di as error "{p 0 0 2}nwdynam requires a network declared via nwset's eventtime() option (an event-type temporal network) - see {help nwset##temporal:nwset}.{p_end}"
		error 198
	}

	// Unit 1 (choice) and Unit 2 (rate, NO intercept) are both
	// implemented - see unw_dynam.do's own header CORRECTION comment
	// for why the rate sub-model turned out to be tractable without
	// real timestamps after all (an ordinal partial likelihood over
	// the sender risk set, same convention as {help nwrem} and this
	// command's own choice sub-model), once no intercept is requested.
	// A WITH-INTERCEPT rate model is a genuinely different, harder
	// continuous-time hazard (confirmed directly from goldfish's own
	// teaching1.Rmd vignette - an intercept-only rate model's own
	// log-likelihood is sensitive to the actual scale of elapsed time)
	// and remains explicitly deferred - there is no way to request an
	// intercept from this command yet, so it is never silently dropped
	// rather than fit. choice_coordination ("expansion batch 11",
	// resolved 2026-09-02, user: "do choice_coordination now") is
	// goldfish's own THIRD DyNAM sub-model - a genuinely different
	// "multinomial-multinomial" joint likelihood over unordered actor
	// pairs (Stadtfeld, Hollway & Block 2017), for undirected
	// ("coordination") tie-formation events, not a flag on choice - see
	// unw_dynam.do's own header comment (above dynam_coord_loglik_
	// grad_multi()) for the formula and its real-goldfish verification.
	if "`submodel'" == "" local submodel "choice"
	if !inlist("`submodel'", "choice", "rate", "choice_coordination") {
		di as error "submodel() must be choice, rate, or choice_coordination - see {help nwdynam}."
		error 198
	}

	// Two-mode (bipartite) networks (added 2026-09-02, "expansion batch
	// 9" - user: "continue"). Previously rejected outright, matching
	// nwrem.ado's own guard - but that guard's own stated reason
	// ("nwset's own twomode/eventtime() composability wasn't available
	// yet") was resolved 2026-08-31 (see the development notes' "Two-mode/
	// temporal architecture initiative"), and unw_dynam.do's own
	// 14-effect Multi engine was never actually checked against a
	// two-mode dyad space until now - it has been, directly against
	// real goldfish (see unw_dynam.do's own header comment on the
	// choice/rate engines for the exact verification and the REAL,
	// goldfish-side per-effect eligibility found - NOT the same as
	// goldfish's own documented `isTwoMode` column, which several
	// direct tests contradicted). goldfish's own architecture assumes
	// a STRICTLY ONE-DIRECTIONAL bipartite dependent network - every
	// event's sender is drawn from node set 1, every receiver from node
	// set 2, never the reverse - matching this package's own
	// contiguous mode-1-then-mode-2 combined actor-index convention
	// (mode 1 = indices 1..n1, mode 2 = n1+1..n, established by
	// nwset's own bipartite ingestion, already reused by nwergm's own
	// bipartite() support). Validated below, once the real event data
	// is available, rather than assumed.
	// choice/rate v1 scope (the development notes scope decision 1,
	// resolved 2026-09-02): directed-only. A two-mode network reports
	// `directed' == false from _nwsyntax (confirmed directly - ERGM's
	// own bipartite convention, where ties have no inherent direction),
	// even though a two-mode DyNAM event network IS inherently
	// directional in the sender/receiver sense goldfish itself assumes
	// (mode 1 always sends, mode 2 always receives) - so the `directed'
	// check is skipped for two-mode networks under choice/rate and real
	// event-level directionality is validated separately below instead,
	// once the actual event data is read.
	// choice_coordination v1 scope (resolved 2026-09-02): the OPPOSITE
	// requirement - a genuinely UNDIRECTED, ONE-MODE network (goldfish's
	// own coordination model is defined over an undirected dependent
	// network; two-mode coordination is a real, disclosed gap - see
	// unw_dynam.do's own header comment on DynamCoordFitMulti for why
	// its own P=p.*p' joint construction needs the actor-1/actor-2
	// candidate pools to be the SAME size, not yet investigated for this
	// package's own two-mode convention).
	if "`submodel'" == "choice_coordination" {
		if "`is2mode'" == "true" {
			di as error "{p 0 0 2}submodel(choice_coordination) does not yet support two-mode (bipartite) networks - see {help nwdynam}.{p_end}"
			error 198
		}
		if "`directed'" == "true" {
			di as error "{p 0 0 2}submodel(choice_coordination) requires an UNDIRECTED network (goldfish's own coordination model is for undirected tie-formation events) - declare `netname' via {opt nwset ..., undirected eventtime()} - see {help nwdynam}.{p_end}"
			error 198
		}
	}
	else if "`directed'" != "true" & "`is2mode'" != "true" {
		di as error "{p 0 0 2}nwdynam requires a directed network for submodel(choice) or submodel(rate) - see {help nwdynam}.{p_end}"
		error 198
	}

	// Effect selection (the development notes' own "effect selection"
	// scope item, resolved 2026-09-02) - each sub-model's own fixed
	// effect list is individually selectable via a boolean flag (or,
	// for the attribute effects added below, a valued option) per
	// effect, matching nwrem.ado's own per-effect-flag convention.
	// `indeg` means something different in each sub-model (the
	// candidate RECEIVER's own in-degree for choice, "alter" type; the
	// candidate SENDER's own in-degree for rate, "ego" type - see
	// {help nwdynam}) - the SAME flag name is reused deliberately
	// (matching goldfish's own indeg() effect, whose type is likewise
	// inferred from context, not a separate option), not two
	// differently-named options for what the reference implementation
	// itself treats as one effect family. Giving NO effect option at
	// all preserves nwdynam's own original behavior exactly (all of
	// the chosen sub-model's STRUCTURAL effects fit together) - flags
	// are for narrowing the set, not a mandatory selection the way
	// {help nwrem}'s own effect options are (nwrem never had an
	// all-effects default to preserve; nwdynam already did, before this
	// option existed, so silently requiring explicit selection now
	// would be a breaking change).
	//
	// same()/diff()/sim()/ego()/alter() (the development notes' "effect
	// expansion" scope item, added 2026-09-02) are VALUED options, not
	// booleans, and - matching goldfish's own effect table exactly -
	// each applies to exactly ONE sub-model: same()/diff()/sim()/alter()
	// choice only (there is no rate-side "same actor as themselves"
	// concept for the dyadic same/diff/sim effects; alter() is
	// dyadic-adjacent - "the candidate's own value" - but still
	// choice-only per goldfish's own table); ego() rate only. `ego()`
	// and `alter()` are the SAME underlying idea (a candidate's own
	// static covariate value, no comparison to anyone) under
	// goldfish's own DIFFERENT names for whichever sub-model's own
	// candidate role it is ("ego" for the rate sub-model's candidate
	// sender; "alter" for the choice sub-model's candidate receiver) -
	// not this package's own naming choice. Unlike the structural
	// flags, none of these five are ever included by the "give
	// nothing, get everything" default (a variable name cannot have a
	// sensible default) - only ever active when a variable is actually
	// named.
	if "`submodel'" == "choice" {
		if "`ego'" != "" {
			di as error "{p 0 0 2}ego() is a rate sub-model effect ({bf:submodel(rate)}), not available with submodel(choice) - see {help nwdynam}.{p_end}"
			error 198
		}
		// outdeg (added 2026-09-02, "expansion batch 5" - user checked
		// this file's own help text against goldfish and found
		// "outdeg: rate sub-model only" WRONG; goldfish's own effect
		// table lists outdeg as valid for BOTH sub-models. `outdeg'
		// means something different in each, exactly like `indeg'
		// already does (candidate RECEIVER's own out-degree here,
		// "alter" type; candidate SENDER's own out-degree under
		// submodel(rate), "ego" type - see unw_dynam.do's own header
		// comment) - the same flag name is reused deliberately, not two
		// differently-named options.
		//
		// trans/cycle/commonsender/commonreceiver/four/nodetrans
		// (added 2026-09-02, "expansion batch 6" - user: "work your way
		// through the list", a systematic audit against goldfish's FULL
		// effect catalog). The first five are genuinely choice-only per
		// goldfish's own table (two-path closure statistics, dyadic);
		// `nodetrans' is - like `indeg'/`outdeg' - valid under BOTH
		// sub-models (same flag name reused deliberately, same
		// convention as those two).
		//
		// egoalterint()/tertius() (added 2026-09-02, "expansion batch
		// 7", continuing the same audit). `egoalterint()' is
		// choice-only (needs both a sender and a candidate role, same
		// structural reason as inertia/recip/same/diff/sim - see
		// unw_dynam.do's own header comment); `tertius()' - like
		// `indeg'/`outdeg'/`nodetrans' - is valid under BOTH sub-models.
		local __nwdynam_effnames "inertia recip indeg same diff sim alter outdeg trans cycle commonsender commonreceiver four nodetrans egoalterint tertius tie"
		local __nwdynam_defaulteffs "inertia recip indeg"
	}
	else if "`submodel'" == "rate" {
		if "`inertia'`recip'`same'`diff'`sim'`alter'`trans'`cycle'`commonsender'`commonreceiver'`four'`egoalterint'`tie'" != "" {
			di as error "{p 0 0 2}inertia/recip/same()/diff()/sim()/alter()/trans/cycle/commonsender/commonreceiver/four/egoalterint()/tie() are choice sub-model effects ({bf:submodel(choice)}), not available with submodel(rate) - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_effnames "indeg outdeg ego nodetrans tertius"
		local __nwdynam_defaulteffs "indeg outdeg"
	}
	else {
		// choice_coordination ("expansion batch 11", 2026-09-02): twelve
		// effects, chosen because each is confirmed BOTH real-
		// goldfish-eligible for choice_coordination (checked directly,
		// not trusted from goldfish's own self-contradicting vignette
		// table - see unw_dynam.do's own header comment on
		// dynam_coord_loglik_grad_multi()) AND already buildable from a
		// reusable state matrix or per-actor covariate in this engine.
		// `recip'/`outdeg'/`cycle'/`commonsender'/`commonreceiver'/
		// `ego' are all confirmed REJECTED by real goldfish for
		// choice_coordination (recip/outdeg structurally make no sense
		// for an undirected tie - no direction to reciprocate, indeg and
		// outdeg coincide). `nodetrans'/`trans' ("batch 12", 2026-09-02),
		// `tertius()'/`four' ("batch 13", 2026-09-02), `egoalterint()'
		// ("batch 14", 2026-09-02, completing the original batch 11
		// eligibility sweep), and `tie()' ("batch 16", 2026-09-02,
		// extending batch 15's own cross-network capability - the
		// simplest possible addition here, since the coordination engine
		// already needs the FULL n x n statistic matrix for every
		// effect, so `S_tie' is literally `tiemat' itself) complete this
		// engine's own currently-wired effect set. `tertiusDiff'/
		// `mixedTrans' are terms of netterms().
		if "`recip'`outdeg'`cycle'`commonsender'`commonreceiver'`ego'" != "" {
			di as error "{p 0 0 2}recip/outdeg/cycle/commonsender/commonreceiver/ego() are not available with submodel(choice_coordination) - only inertia/indeg/same()/diff()/sim()/alter()/nodetrans/trans/tertius()/four/egoalterint()/tie() are wired for this sub-model so far (tertiusdiff and mixedtrans are terms of netterms()) - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_effnames "inertia indeg same diff sim alter nodetrans trans tertius four egoalterint tie"
		local __nwdynam_defaulteffs "inertia indeg"
	}

	// Two-mode effect eligibility (added 2026-09-02, "expansion batch
	// 9"). goldfish's own DOCUMENTED `isTwoMode` column (per its own
	// doc/goldfishEffects.Rmd) turned out NOT to reliably predict which
	// effects actually work on a two-mode network - direct testing
	// against real goldfish found several effects the table marks
	// eligible (`recip`, `outdeg` choice-side, `commonReceiver`,
	// `indeg` rate-side) hard-REJECTED by goldfish's own engine at
	// runtime, each with a genuine structural reason: goldfish's own
	// two-mode DyNAM architecture is STRICTLY one-directional (mode 1
	// only ever sends, mode 2 only ever receives), so any effect that
	// would need mode 1's own in-ties or mode 2's own out-ties is
	// undefined (there are none, by construction) - `recip` (does the
	// candidate have a tie BACK to the sender - impossible, mode 2
	// never sends), `outdeg` choice-side (candidate's own out-degree -
	// impossible, mode 2 never sends), `commonReceiver`/`indeg`
	// rate-side both need the SAME missing role. A few effects goldfish
	// DOES accept for two-mode (`commonSender`, `four`) gave
	// non-degenerate fits that naive risk-set-only reasoning could not
	// explain, matching goldfish's own doc warning ("we omit the
	// difference in the computation of the statistics when isTwoMode
	// is used") - some effects genuinely redefine their own formula for
	// two-mode. This was confirmed to run DEEPER than expected: direct
	// comparison against goldfish's own reconstructed statistics
	// (`preprocessingOnly=TRUE`, same technique as the `four()`/
	// `tertius()` fixes elsewhere in this file) found `same()` and
	// `ego()` ALSO silently wrong under the naive combined-vector
	// approach originally attempted here - goldfish's own internal
	// two-mode array applies its usual "zero the self-tie" convention
	// by raw ROW-INDEX-EQUALS-COLUMN-INDEX matching even though row and
	// column span two GENUINELY DIFFERENT node sets in a two-mode
	// network (so e.g. person-index-1 vs org-index-1 gets zeroed as if
	// they were "the same actor," despite being different actors in
	// different modes) - a real goldfish-internal representation detail
	// that could not be confidently reverse-engineered in the time
	// available. Given this, EVERY attribute effect (`same`/`diff`/
	// `sim`/`alter`/`ego`/`egoalterint`/`tertius`) is rejected for
	// two-mode, not just the structurally-impossible ones - disclosed
	// as a genuine unresolved gap, not silently shipped with an
	// unverified formula. Verified WORKING and implemented, with an
	// EXACT numerical match to real goldfish: `inertia`, `indeg`
	// (choice, alter type), `outdeg` (rate, ego type) - see
	// unw_dynam.do's own header comment for the numbers.
	// `nodeTrans`/`trans`/`cycle` are excluded per goldfish's own
	// `isTwoMode=×` table entries (not independently re-verified, since
	// every OTHER table entry that WAS re-verified turned out
	// unreliable in at least one direction - disclosed, not assumed
	// safe either way). Windowed and weighted effects are not yet
	// verified for two-mode at all - rejected outright rather than
	// silently combined.
	if "`is2mode'" == "true" {
		if "`submodel'" == "choice" {
			local __nwdynam_defaulteffs "inertia indeg"
			if "`recip'" != "" {
				di as error "{p 0 0 2}recip is not supported for two-mode (bipartite) networks - goldfish itself rejects it (mode 2 never sends, so a candidate can never have a tie back to the sender) - see {help nwdynam}.{p_end}"
				error 198
			}
			if "`outdeg'" != "" {
				di as error "{p 0 0 2}outdeg is not supported for two-mode (bipartite) networks under submodel(choice) - goldfish itself rejects it (mode 2 never sends, so a candidate has no out-degree) - see {help nwdynam}.{p_end}"
				error 198
			}
			foreach __nwdynam_e in trans cycle commonsender commonreceiver nodetrans same diff sim {
				if "``__nwdynam_e''" != "" {
					di as error "{p 0 0 2}`__nwdynam_e' is not available for two-mode (bipartite) networks - goldfish itself rejects it or computes it with the wrong dimensions there - see {help nwdynam}.{p_end}"
					error 198
				}
			}
		}
		else {
			local __nwdynam_defaulteffs "outdeg"
			if "`indeg'" != "" {
				di as error "{p 0 0 2}indeg is not supported for two-mode (bipartite) networks under submodel(rate) - goldfish itself rejects it (mode 1 never receives, so a candidate sender has no in-degree) - see {help nwdynam}.{p_end}"
				error 198
			}
			foreach __nwdynam_e in nodetrans tertius {
				if "``__nwdynam_e''" != "" {
					di as error "{p 0 0 2}`__nwdynam_e' is not available for two-mode (bipartite) networks under submodel(rate) - goldfish itself rejects it - see {help nwdynam}.{p_end}"
					error 198
				}
			}
		}
		if `inertiawindow' != -1 | `recipwindow' != -1 | `indegwindow' != -1 | `outdegwindow' != -1 {
			di as error "{p 0 0 2}windowed effects are not yet verified for two-mode (bipartite) networks - see {help nwdynam}.{p_end}"
			error 198
		}
		if "`weightedinertia'`weightedrecip'`weightedindeg'`weightedoutdeg'" != "" {
			di as error "{p 0 0 2}weighted effects are not yet verified for two-mode (bipartite) networks - see {help nwdynam}.{p_end}"
			error 198
		}
		if "`change'" != "" & "`egoalterint'" != "" {
			di as error "{p 0 0 2}change() together with egoalterint() is not available for two-mode (bipartite) networks: goldfish does not give a usable result for this combination, so it stays refused.{p_end}"
			error 198
		}
		if "`tie'`opportunities'" != "" {
			di as error "{p 0 0 2}tie() and opportunities() are not yet verified for two-mode (bipartite) networks - see {help nwdynam}.{p_end}"
			error 198
		}
	}

	// inertiawindow()/recipwindow() (the development notes' "effect
	// expansion" scope item, batch 3, added 2026-09-02 - user request
	// "implement window effects," REDESIGNED twice after user feedback:
	// first to two independent plain options (a nested suboption
	// syntax, e.g. inertia(window(#)), "does not feel like Stata or
	// the rest of nwcommands"), then a second time so that GIVING the
	// window option is itself sufficient to activate its own effect -
	// "handle it how we handled gwesp in nwergm... recipwindow would
	// not require recip." Matching {help nwergm}'s own `gwesp(real)`
	// convention exactly: there is no separate "include gwesp" flag,
	// giving `gwesp(0.5)` both selects the term AND sets its own
	// parameter in one step - `inertiawindow(#)`/`recipwindow(#)` now
	// work the same way for `inertia`/`recip`. `-1` sentinel = "not
	// given," matching {opt seed()}'s own established convention in
	// this same file.
	// netterms(): effects on other networks, mixed closures, tertiusDiff (choice and choice_coordination). tie() of an event network
	// (a network that changes over time) is such a term. The terms are parsed here, the networks resolved further down.
	tempname __nwdynam_xt __nwdynam_XL __nwdynam_XA
	local __nwdynam_nxt = 0
	capture mata: mata drop __nwdynam_xblob
	if `"`tie'"' != "" & "`submodel'" != "rate" & "`is2mode'" != "true" {
		_nwsyntax `tie', other(__nwdynam_tienet_) max(1)
		if "`__nwdynam_tienet_istemporal'" == "true" {
			local netterms `"tie(`tie') `netterms'"'
			local tie ""
		}
	}
	if `"`netterms'"' != "" {
		if "`is2mode'" == "true" {
			di as error "{p 0 0 2}netterms() is not yet available for two-mode (bipartite) networks - see {help nwdynam}.{p_end}"
			error 198
		}
		mata: `__nwdynam_xt' = dynam_parse_terms(`"`netterms'"', "`submodel'", "`netname'")
		mata: st_local("__nwdynam_nxt", strofreal(rows(`__nwdynam_xt')))
	}

	local __nwdynam_windowinertia = 1e300
	local __nwdynam_windowrecip = 1e300
	if `inertiawindow' != -1 {
		if "`submodel'" == "rate" {
			di as error "{p 0 0 2}inertiawindow() only applies to the inertia effect of the choice and choice_coordination sub-models - see {help nwdynam}.{p_end}"
			error 198
		}
		if `inertiawindow' <= 0 {
			di as error "{p 0 0 2}inertiawindow() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
			error 198
		}
		local __nwdynam_windowinertia = `inertiawindow'
		local inertia "inertia"
	}
	if `recipwindow' != -1 {
		if "`submodel'" != "choice" {
			di as error "{p 0 0 2}recipwindow() only applies to the choice sub-model's recip effect - see {help nwdynam}.{p_end}"
			error 198
		}
		if `recipwindow' <= 0 {
			di as error "{p 0 0 2}recipwindow() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
			error 198
		}
		local __nwdynam_windowrecip = `recipwindow'
		local recip "recip"
	}

	// indegwindow()/outdegwindow() (added 2026-09-02, user follow-up
	// "can we make [window] also available for the rate part of the
	// model?", immediately after inertiawindow()/recipwindow()'s own
	// gwesp-style self-activation redesign): the SAME self-activating
	// convention, for the rate sub-model's own indeg/outdeg effects.
	// goldfish's own doc/goldfishEffects.Rmd effect table confirms
	// `window` is a real argument on indeg/outdeg and that both apply to
	// the rate sub-model (see unw_dynam.do's own
	// dynam_rate_loglik_grad_multi() header comment for the exact
	// formula and its real-goldfish verification). Rate-only, matching
	// inertiawindow()/recipwindow()'s own choice-only restriction (the
	// SAME flag name `indeg` already means a different thing per
	// sub-model - see the effect-selection comment above - so
	// indegwindow() is unambiguous once restricted to submodel(rate)).
	local __nwdynam_windowindeg = 1e300
	local __nwdynam_windowoutdeg = 1e300
	if `indegwindow' != -1 {
		if "`submodel'" == "choice" {
			di as error "{p 0 0 2}indegwindow() only applies to the indeg effect of the rate and choice_coordination sub-models - see {help nwdynam}.{p_end}"
			error 198
		}
		if `indegwindow' <= 0 {
			di as error "{p 0 0 2}indegwindow() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
			error 198
		}
		local __nwdynam_windowindeg = `indegwindow'
		local indeg "indeg"
	}
	if `outdegwindow' != -1 {
		if "`submodel'" != "rate" {
			di as error "{p 0 0 2}outdegwindow() only applies to the rate sub-model's outdeg effect - see {help nwdynam}.{p_end}"
			error 198
		}
		if `outdegwindow' <= 0 {
			di as error "{p 0 0 2}outdegwindow() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
			error 198
		}
		local __nwdynam_windowoutdeg = `outdegwindow'
		local outdeg "outdeg"
	}

	// tertiuswindow(#): window of the tertius() effect (all sub-models): only the ties of the last # time units count
	local __nwdynam_windowtertius = 1e300
	if `tertiuswindow' != -1 {
		if "`tertius'" == "" {
			di as error "{p 0 0 2}tertiuswindow() needs tertius(): it is the window of that effect.{p_end}"
			error 198
		}
		if `tertiuswindow' <= 0 {
			di as error "{p 0 0 2}tertiuswindow() must be a positive number (in the same time units as the network's own eventtime()).{p_end}"
			error 198
		}
		local __nwdynam_windowtertius = `tertiuswindow'
	}

	// weightedinertia/weightedrecip/weightedindeg/weightedoutdeg (added
	// 2026-09-02, "expansion batch 8" - user: "work your way through
	// the list"): goldfish's own `weighted=TRUE' argument on
	// inertia/recip/indeg/outdeg - the effect reads the value of the
	// tie (the sum of the values of its prior events, 1 per event when
	// the network has no event values) instead of binary presence.
	// With a window the value within the window, as goldfish. Self-activating, matching inertiawindow()/recipwindow()'s
	// own gwesp-style convention exactly - no separate boolean flag
	// needed. `weightedindeg'/`weightedoutdeg' are valid under BOTH
	// sub-models (goldfish's own table confirms this, unlike
	// indegwindow()/outdegwindow() which remain rate-only - a real,
	// disclosed choice-side windowing gap, not a weighted one).
	local __nwdynam_weightedinertia = 0
	local __nwdynam_weightedrecip = 0
	local __nwdynam_weightedindeg = 0
	local __nwdynam_weightedoutdeg = 0
	if "`weightedinertia'" != "" {
		if "`submodel'" == "rate" {
			di as error "{p 0 0 2}weightedinertia only applies to the inertia effect of the choice and choice_coordination sub-models - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_weightedinertia = 1
		local inertia "inertia"
	}
	if "`weightedrecip'" != "" {
		if "`submodel'" != "choice" {
			di as error "{p 0 0 2}weightedrecip only applies to the choice sub-model's recip effect - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_weightedrecip = 1
		local recip "recip"
	}
	if "`weightedindeg'" != "" {
		local __nwdynam_weightedindeg = 1
		local indeg "indeg"
	}
	if "`weightedoutdeg'" != "" {
		if "`submodel'" == "choice_coordination" {
			di as error "{p 0 0 2}weightedoutdeg is not available with submodel(choice_coordination) - outdeg itself is not either (indeg and outdeg coincide for an undirected network) - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_weightedoutdeg = 1
		local outdeg "outdeg"
	}

	// intercept ("expansion batch 17", 2026-09-02, continuing the
	// goldfish-comparison list) - the genuinely continuous-time
	// competing-risks hazard variant this file's own header comment
	// always described as deferred. Unlike the NO-intercept rate model
	// (which reduces to the same ordinal partial likelihood as
	// everything else in this package, using only event ORDER, never
	// real elapsed time), giving `intercept' makes the model genuinely
	// sensitive to the ACTUAL SCALE of elapsed real time between events
	// (`netname''s own eventtime() values) - matching goldfish's own
	// `dep ~ 1 + ...' convention (a literal `1' as the first formula
	// term requests an intercept; verified directly against goldfish's
	// own `parseIntercept()' source). submodel(rate) only - matching
	// goldfish's own DyNAM architecture, where only the rate sub-model
	// has ever had an intercept concept (the choice sub-model's own
	// conditional logit has no intercept term at all, by construction).
	// v1 scope, extended "expansion batch 18" (2026-09-02): two-mode and
	// `weighted=TRUE' are now verified together with `intercept' (each
	// independently, against real goldfish on a toy network with real
	// timestamps - dev/dynam_unit21_rateintercept_twomode_crosscheck.R,
	// dynam_unit21c_..._weighted_...R) and no longer rejected - the SAME
	// statistic swaps the no-intercept engine already uses (a real-time
	// window read, a cumulative weighted count, mode-2 actors masked
	// out of the hazard sum) turned out to compose correctly with the
	// hazard-integral aggregation exactly as written, confirmed by
	// exact agreement with real goldfish, not merely assumed.
	//
	// `indegwindow()'/`outdegwindow()' remain REJECTED with `intercept'
	// - this one is a genuine, disclosed architectural gap, not a
	// pending verification: evaluating the fitted log-likelihood
	// directly at real goldfish's own reported MLE (not re-optimizing)
	// found a real, substantial mismatch (~0.3 in logLik, far past
	// floating-point noise), traced to goldfish's own `window()'
	// mechanism (`createWindowedEvents()' in its R source) inserting
	// SYNTHETIC "dissolve" events into the timeline at exactly
	// `event_time + window' for every windowed tie, forcing its
	// piecewise-constant hazard machinery to recompute right at each
	// expiry. This engine only recomputes statistics at REAL dependent
	// events, so a windowed contact silently keeps counting toward
	// `indeg'/`outdeg' for the rest of the current inter-event interval
	// even after it should have aged out mid-interval - wrong only when
	// intercept's own hazard-INTEGRAL aggregation is in play (the
	// no-intercept engine's ordinal partial likelihood never integrates
	// over real time at all, which is why window already works
	// correctly there). Fixing this needs the same synthetic-event
	// injection goldfish itself uses - a real, disclosed follow-on, not
	// attempted this batch.
	if "`intercept'" != "" {
		if "`submodel'" != "rate" {
			di as error "intercept only applies to submodel(rate) - see {help nwdynam}."
			error 198
		}
	}

	local __nwdynam_anygiven = 0
	foreach __nwdynam_e of local __nwdynam_effnames {
		if "``__nwdynam_e''" != "" local __nwdynam_anygiven = 1
	}
	if `__nwdynam_nxt' > 0 local __nwdynam_anygiven = 1
	// ignorerep(): effects whose statistic is set to 0 for the pairs that already have a tie in the dependent network (goldfish: ignoreRep = TRUE)
	if "`ignorerep'" != "" {
		if "`submodel'" == "rate" {
			di as error "{p 0 0 2}ignorerep() is available for submodel(choice) and submodel(choice_coordination); the rate sub-model is not available for ignorerep().{p_end}"
			error 198
		}
		foreach __nwdynam_w of local ignorerep {
			local __nwdynam_w = lower("`__nwdynam_w'")
			if inlist("`__nwdynam_w'", "same", "diff", "sim", "alter", "ego", "egoalterint") {
				di as error "{p 0 0 2}ignorerep(): {bf:`__nwdynam_w'} has no ignoreRep in goldfish (only structural effects have).{p_end}"
				error 198
			}
			local __nwdynam_ok : list posof "`__nwdynam_w'" in __nwdynam_effnames
			if `__nwdynam_ok' == 0 {
				di as error "{p 0 0 2}ignorerep(): {bf:`__nwdynam_w'} is not an effect of submodel(`submodel'); give the effects of the model (the terms of netterms() have their own ignorerep option).{p_end}"
				error 198
			}
			if "``__nwdynam_w''" == "" & !(`__nwdynam_anygiven' == 0 & `: list posof "`__nwdynam_w'" in __nwdynam_defaulteffs' > 0) {
				di as error "{p 0 0 2}ignorerep(): {bf:`__nwdynam_w'} is not among the effects of this model.{p_end}"
				error 198
			}
		}
		local ignorerep = lower("`ignorerep'")
	}

	local __nwdynam_active ""
	local __nwdynam_sep ""
	local __nwdynam_collabs ""
	local __nwdynam_nactive = 0
	foreach __nwdynam_e of local __nwdynam_effnames {
		local __nwdynam_isdefault : list posof "`__nwdynam_e'" in __nwdynam_defaulteffs
		if "``__nwdynam_e''" != "" | (`__nwdynam_anygiven' == 0 & `__nwdynam_isdefault' > 0) {
			local __nwdynam_ig : list posof "`__nwdynam_e'" in ignorerep
			local __nwdynam_active "`__nwdynam_active'`__nwdynam_sep'`=cond(`__nwdynam_ig' > 0, 2, 1)'"
			local __nwdynam_collabs "`__nwdynam_collabs' `__nwdynam_e'"
			local __nwdynam_nactive = `__nwdynam_nactive' + 1
		}
		else {
			local __nwdynam_active "`__nwdynam_active'`__nwdynam_sep'0"
		}
		local __nwdynam_sep ","
	}

	// same()/diff()/sim()/ego()/alter()/tertius()/egoalterint() covariates are read from the current dataset. When it
	// carries the node labels (variable _nwnode, as after nwload, xvars), every actor is matched to its row by label,
	// whatever the sort order; otherwise the data must have one row per actor, in actor order. Missing values are
	// refused.
	tempname __nwdynam_same __nwdynam_diff __nwdynam_sim __nwdynam_ego __nwdynam_alter __nwdynam_tertius __nwdynam_rows
	tempname __nwdynam_egoalterint1 __nwdynam_egoalterint2
	tempname __nwdynam_pinit
	local __nwdynam_xvars ""
	if `__nwdynam_nxt' > 0 mata: st_local("__nwdynam_xvars", invtokens(dynam_xt_vars(`__nwdynam_xt')))
	foreach __nwdynam_xv of local __nwdynam_xvars {
		local __nwdynam_ok : list posof "`__nwdynam_xv'" in __nwdynam_cvars
		if `__nwdynam_ok' > 0 {
			di as error "{p 0 0 2}change(): the covariate {bf:`__nwdynam_xv'} is used by tertiusdiff() in netterms(); a change of a covariate is not available for tertiusdiff().{p_end}"
			error 198
		}
	}
	if "`same'`diff'`sim'`ego'`alter'`tertius'`egoalterint'`present'`__nwdynam_xvars'" != "" {
		unw_defs
		capture confirm string variable `nw_nodename'
		if !_rc {
			mata: `__nwdynam_rows' = dynam_node_rows(`netobj'->get_nodenames(), st_sdata(., "`nw_nodename'"))
			mata: st_local("__nwdynam_nf", strofreal(sum(`__nwdynam_rows' :== 0)))
			mata: st_local("__nwdynam_nd", strofreal(sum(`__nwdynam_rows' :< 0)))
			if `__nwdynam_nf' > 0 {
				mata: st_local("__nwdynam_lab", select(`netobj'->get_nodenames(), (`__nwdynam_rows' :== 0)')[1])
				capture mata: mata drop `__nwdynam_rows'
				di as error "{p 0 0 2}`__nwdynam_nf' actor(s) of `netname' (first: {bf:`__nwdynam_lab'}) have no row in the current dataset (variable `nw_nodename'); covariates need one row per actor (see {help nwload}'s {bf:xvars} option).{p_end}"
				error 198
			}
			if `__nwdynam_nd' > 0 {
				mata: st_local("__nwdynam_lab", select(`netobj'->get_nodenames(), (`__nwdynam_rows' :< 0)')[1])
				capture mata: mata drop `__nwdynam_rows'
				di as error "{p 0 0 2}`__nwdynam_nd' actor(s) of `netname' (first: {bf:`__nwdynam_lab'}) have more than one row in the current dataset (variable `nw_nodename'); covariates need exactly one row per actor.{p_end}"
				error 198
			}
		}
		else {
			if _N != `nodes' {
				di as error "{p 0 0 2}the current dataset has `=_N' observations but `netname' has `nodes' actors and no node labels (variable `nw_nodename'); covariates need one row per actor, in actor order (see {help nwload}'s {bf:xvars} option).{p_end}"
				error 198
			}
			mata: `__nwdynam_rows' = (1::`nodes')
		}
	}
	local __nwdynam_ea1 : word 1 of `egoalterint'
	local __nwdynam_ea2 : word 2 of `egoalterint'
	foreach __nwdynam_cv in same diff sim ego alter tertius egoalterint1 egoalterint2 {
		if "`__nwdynam_cv'" == "egoalterint1" local __nwdynam_var "`__nwdynam_ea1'"
		else if "`__nwdynam_cv'" == "egoalterint2" local __nwdynam_var "`__nwdynam_ea2'"
		else local __nwdynam_var "``__nwdynam_cv''"
		if "`__nwdynam_var'" != "" {
			mata: `__nwdynam_`__nwdynam_cv'' = st_data(`__nwdynam_rows', "`__nwdynam_var'")'
			mata: st_local("__nwdynam_nmiss", strofreal(missing(`__nwdynam_`__nwdynam_cv'')))
			if `__nwdynam_nmiss' > 0 {
				capture mata: mata drop `__nwdynam_rows' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
				di as error "{p 0 0 2}covariate `__nwdynam_var': `__nwdynam_nmiss' actor(s) of `netname' have a missing value; the covariate must be defined for every actor.{p_end}"
				error 416
			}
		}
		else {
			mata: `__nwdynam_`__nwdynam_cv'' = J(1, `nodes', 0)
		}
	}
	if "`present'" != "" {
		mata: `__nwdynam_pinit' = st_data(`__nwdynam_rows', "`present'")'
		mata: st_local("__nwdynam_nmiss", strofreal(missing(`__nwdynam_pinit')))
		mata: st_local("__nwdynam_nbad", strofreal(sum((`__nwdynam_pinit' :!= 0) :& (`__nwdynam_pinit' :!= 1))))
		if `__nwdynam_nmiss' > 0 | `__nwdynam_nbad' > 0 {
			capture mata: mata drop `__nwdynam_rows' `__nwdynam_pinit' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
			di as error "{p 0 0 2}present(): variable `present' must be 1 (present at the start) or 0 (absent at the start) for every actor of `netname', with no missing values.{p_end}"
			error 198
		}
	}
	else mata: `__nwdynam_pinit' = J(1, 0, 0)
	mata: `__nwdynam_XA' = asarray_create("string", 1)
	foreach __nwdynam_xv of local __nwdynam_xvars {
		capture confirm numeric variable `__nwdynam_xv'
		if _rc {
			capture mata: mata drop `__nwdynam_rows'
			di as error "{p 0 0 2}netterms(): {bf:`__nwdynam_xv'} is not a numeric variable of the current dataset.{p_end}"
			error 109
		}
		mata: asarray(`__nwdynam_XA', "`__nwdynam_xv'", st_data(`__nwdynam_rows', "`__nwdynam_xv'")')
		mata: st_local("__nwdynam_nmiss", strofreal(missing(asarray(`__nwdynam_XA', "`__nwdynam_xv'"))))
		if `__nwdynam_nmiss' > 0 {
			capture mata: mata drop `__nwdynam_rows'
			di as error "{p 0 0 2}covariate `__nwdynam_xv': `__nwdynam_nmiss' actor(s) of `netname' have a missing value; the covariate must be defined for every actor.{p_end}"
			error 416
		}
	}
	capture mata: mata drop `__nwdynam_rows'

	tempname __nwdynam_chg
	mata: `__nwdynam_chg' = J(0, 4, 0)
	forvalues __nwdynam_i = 1/`__nwdynam_nchgpairs' {
		local __nwdynam_cvar : word `__nwdynam_i' of `__nwdynam_cvars'
		local __nwdynam_cnet : word `__nwdynam_i' of `__nwdynam_cnets'
		_nwsyntax `__nwdynam_cnet', other(__nwdynam_chgnet_) max(1)
		if "`__nwdynam_chgnet_istemporal'" != "true" | "`__nwdynam_chgnet_temporaltype'" != "event" {
			capture mata: mata drop `__nwdynam_chg' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
			di as error "{p 0 0 2}change(): `__nwdynam_cnet' is not an event network; declare the changes with {bf:nwset id id value, eventtime(time) selfloop name(`__nwdynam_cnet')}.{p_end}"
			error 198
		}
		// the covariate rows that use this variable: 1 same 2 diff 3 sim 4 ego/alter 5 tertius 6/7 egoalterint
		local __nwdynam_slots ""
		if "`same'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 1"
		if "`diff'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 2"
		if "`sim'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 3"
		if "`ego'`alter'" != "" & ("`ego'" == "`__nwdynam_cvar'" | "`alter'" == "`__nwdynam_cvar'") local __nwdynam_slots "`__nwdynam_slots' 4"
		if "`tertius'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 5"
		if "`__nwdynam_ea1'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 6"
		if "`__nwdynam_ea2'" == "`__nwdynam_cvar'" local __nwdynam_slots "`__nwdynam_slots' 7"
		local __nwdynam_firstslot : word 1 of `__nwdynam_slots'
		local __nwdynam_initvec : word `__nwdynam_firstslot' of same diff sim cov4 tertius egoalterint1 egoalterint2
		if `__nwdynam_firstslot' == 4 {
			if "`ego'" != "" mata: __nwdynam_iv = `__nwdynam_ego'
			else mata: __nwdynam_iv = `__nwdynam_alter'
		}
		else mata: __nwdynam_iv = `__nwdynam_`__nwdynam_initvec''
		mata: __nwdynam_cr = *(`__nwdynam_chgnet_netobj'->get_eventlist())
		capture noisily mata: `__nwdynam_chg' = `__nwdynam_chg' \ dynam_change_rows(__nwdynam_cr, `__nwdynam_chgnet_netobj'->get_nodenames(), `netobj'->get_nodenames(), __nwdynam_iv, ("`increment'" != ""), strtoreal(tokens("`__nwdynam_slots'")), "`__nwdynam_cnet'", "`__nwdynam_cvar'")
		local __nwdynam_crc = _rc
		capture mata: mata drop __nwdynam_iv __nwdynam_cr
		if `__nwdynam_crc' {
			capture mata: mata drop `__nwdynam_chg' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
			exit `__nwdynam_crc'
		}
	}
	if `__nwdynam_nchgpairs' > 1 {
		// the order inside a time point: event list order of the pairs, as given; sort by time, stable
		mata: `__nwdynam_chg' = `__nwdynam_chg'[order((`__nwdynam_chg'[., 3], (1::rows(`__nwdynam_chg'))), (1, 2)), .]
	}
	// composition change: rows with code 8 (join/leave events) and 9 (presence at the start) after the covariate changes
	local __nwdynam_ncomp = 0
	tempname __nwdynam_comprows
	mata: `__nwdynam_comprows' = J(0, 4, 0)
	if `"`compchange'"' != "" {
		if `: word count `compchange'' != 1 {
			capture mata: mata drop `__nwdynam_comprows' `__nwdynam_pinit' `__nwdynam_chg' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
			di as error "{p 0 0 2}compchange() takes the name of one event network with the joins (value 1) and leaves (value 0).{p_end}"
			error 198
		}
		_nwsyntax `compchange', other(__nwdynam_compnet_) max(1)
		if "`__nwdynam_compnet_istemporal'" != "true" | "`__nwdynam_compnet_temporaltype'" != "event" {
			capture mata: mata drop `__nwdynam_comprows' `__nwdynam_pinit' `__nwdynam_chg' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
			di as error "{p 0 0 2}compchange(): `compchange' is not an event network; declare the joins and leaves with {bf:nwset id id value, eventtime(time) selfloop name(`compchange')}.{p_end}"
			error 198
		}
		mata: __nwdynam_cr = *(`__nwdynam_compnet_netobj'->get_eventlist())
	}
	else mata: __nwdynam_cr = J(0, 4, 0)
	local __nwdynam_cnames `"J(1, 0, "")"'
	if `"`compchange'"' != "" local __nwdynam_cnames "`__nwdynam_compnet_netobj'->get_nodenames()"
	capture noisily mata: `__nwdynam_comprows' = dynam_comp_rows(__nwdynam_cr, `__nwdynam_cnames', `netobj'->get_nodenames(), `__nwdynam_pinit', "`compchange'")
	local __nwdynam_crc = _rc
	capture mata: mata drop __nwdynam_cr `__nwdynam_pinit'
	if `__nwdynam_crc' {
		capture mata: mata drop `__nwdynam_comprows' `__nwdynam_chg' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2'
		exit `__nwdynam_crc'
	}
	mata: st_local("__nwdynam_ncomp", strofreal(sum(`__nwdynam_comprows'[., 1] :== 8)))

	if `"`seed'"' != "" {
		set seed `seed'
	}

	if `__nwdynam_nxt' > 0 {
		mata: st_local("__nwdynam_xnets", invtokens(dynam_xt_nets(`__nwdynam_xt', "`netname'")))
		mata: `__nwdynam_XL' = asarray_create("string", 1)
		foreach __nwdynam_xn of local __nwdynam_xnets {
			_nwsyntax `__nwdynam_xn', other(__nwdynam_xnet_) max(1)
			if "`__nwdynam_xnet_istemporal'" == "true" & "`__nwdynam_xnet_temporaltype'" != "event" {
				di as error "{p 0 0 2}netterms(): `__nwdynam_xn' is a temporal network that is not an event network; use an ordinary network or one declared with {bf:nwset ..., eventtime()}.{p_end}"
				error 198
			}
			if "`__nwdynam_xnet_istemporal'" == "true" {
				mata: asarray(`__nwdynam_XL', "`__nwdynam_xn'", dynam_xlayer("true", "`__nwdynam_xnet_directed'", "`__nwdynam_xn'", `netobj'->get_nodenames(), *(`__nwdynam_xnet_netobj'->get_eventlist()), `__nwdynam_xnet_netobj'->get_nodenames(), J(0, 0, 0)))
			}
			else {
				mata: asarray(`__nwdynam_XL', "`__nwdynam_xn'", dynam_xlayer("false", "`__nwdynam_xnet_directed'", "`__nwdynam_xn'", `netobj'->get_nodenames(), J(0, 0, 0), J(1, 0, ""), *(`__nwdynam_xnet_netobj'->get_matrix())))
			}
		}
		mata: __nwdynam_xblob = dynam_xblob(`__nwdynam_xt', tokens("`__nwdynam_xnets'"), `__nwdynam_XL', "`netname'", tokens("`__nwdynam_xvars'"), `__nwdynam_XA', `nodes')
		mata: st_local("__nwdynam_xlabs", invtokens(dynam_xt_labels(`__nwdynam_xt', "`netname'")))
		local __nwdynam_collabs "`__nwdynam_collabs' `__nwdynam_xlabs'"
	}
	capture mata: mata drop `__nwdynam_xt' `__nwdynam_XL' `__nwdynam_XA'

	tempname __nwdynam_events
	mata: `__nwdynam_events' = *(`netobj'->get_eventlist())
	mata: st_numscalar("__nwdynam_nevents", rows(`__nwdynam_events'))
	local nevents = __nwdynam_nevents
	scalar drop __nwdynam_nevents
	// Events with an identical time are simultaneous (as in goldfish): statistics for all events of a time point are
	// computed from the state before that time point. Option sequential: old behaviour, tied events one after another.
	mata: st_numscalar("__nwdynam_ntimes", rows(uniqrows(`__nwdynam_events'[.,3])))
	local __nwdynam_ntimes = __nwdynam_ntimes
	scalar drop __nwdynam_ntimes

	// events whose sender or receiver is not an actor of the network (a missing or unmatched id; actor index 0) or
	// that have no time cannot be placed in the event sequence: refuse them
	mata: st_local("__nwdynam_bad0", strofreal(sum((`__nwdynam_events'[., 1] :< 1) :| (`__nwdynam_events'[., 2] :< 1) :| missing(`__nwdynam_events'[., 1]) :| missing(`__nwdynam_events'[., 2]))))
	mata: st_local("__nwdynam_badt", strofreal(missing(`__nwdynam_events'[., 3])))
	mata: st_local("__nwdynam_badv", strofreal(cols(`__nwdynam_events') >= 4 ? missing(`__nwdynam_events'[., cols(`__nwdynam_events') >= 4 ? 4 : 3]) : 0))
	if `__nwdynam_bad0' > 0 | `__nwdynam_badt' > 0 | `__nwdynam_badv' > 0 {
		_nwdynam_dropmata `__nwdynam_events' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2' `__nwdynam_chg' `__nwdynam_comprows'
		if `__nwdynam_bad0' > 0 {
			di as error "{p 0 0 2}`__nwdynam_bad0' event(s) of `netname' have a sender or receiver that is not an actor of the network (a missing actor id, or an id that does not match a node label); drop these events before {cmd:nwset}.{p_end}"
		}
		if `__nwdynam_badt' > 0 {
			di as error "{p 0 0 2}`__nwdynam_badt' event(s) of `netname' have a missing event time; drop these events before {cmd:nwset}.{p_end}"
		}
		if `__nwdynam_badv' > 0 {
			di as error "{p 0 0 2}`__nwdynam_badv' event(s) of `netname' have a missing value; drop these events before {cmd:nwset}.{p_end}"
		}
		error 416
	}

	if `nevents' < 2 {
		di as error "nwdynam needs at least 2 events to fit a model; `netname' has `nevents'."
		error 2001
	}

	// opportunities() ("expansion batch 10", added 2026-09-02, user:
	// "work your way through it" - continuing the goldfish-comparison
	// list) - goldfish's own `opportunitiesList' estimationInit
	// argument, choice sub-model only (goldfish's own documentation:
	// "ONLY for choice models"). Takes exactly TWO variables from the
	// CURRENT Stata dataset: an event sequence number (1..`nevents',
	// matching `netname''s own chronological event order) and an actor
	// ID (1..`nodes'), ONE ROW PER (event, available-actor) PAIR - a
	// genuinely different dataset SHAPE from same()/diff()/sim()/
	// alter()/ego()/tertius()/egoalterint() (one row per ACTOR) or the
	// event-level dataset `netname' itself was declared from, so
	// opportunities() cannot be combined with those seven options in
	// the same call (a real, disclosed v1 limitation - the current
	// dataset cannot simultaneously have both shapes). Mechanics
	// (numeric node indices, not labels; the sender's own
	// self-exclusion still applies regardless of whether it is itself
	// listed) derived from a real toy example against goldfish FIRST -
	// see unw_dynam.do's own header comment and
	// dev/dynam_unit13_opportunities_crosscheck.R/.do for the
	// real-goldfish verification.
	tempname __nwdynam_oppmat
	if "`opportunities'" != "" {
		if "`submodel'" != "choice" {
			di as error "{p 0 0 2}opportunities() only applies to the choice sub-model (goldfish's own opportunitiesList is choice-only) - see {help nwdynam}.{p_end}"
			error 198
		}
		local __nwdynam_oppev : word 1 of `opportunities'
		local __nwdynam_oppact : word 2 of `opportunities'
		local __nwdynam_fpre ""
		if "`oppframe'" != "" local __nwdynam_fpre "frame `oppframe':"
		quietly `__nwdynam_fpre' count if missing(`__nwdynam_oppev') | missing(`__nwdynam_oppact')
		if r(N) > 0 {
			di as error "{p 0 0 2}opportunities(): `r(N)' row(s) have a missing value in `__nwdynam_oppev' or `__nwdynam_oppact'; drop them first.{p_end}"
			error 416
		}
		capture `__nwdynam_fpre' assert `__nwdynam_oppev' >= 1 & `__nwdynam_oppev' <= `nevents' if !missing(`__nwdynam_oppev')
		if _rc {
			di as error "{p 0 0 2}`__nwdynam_oppev' must be between 1 and `nevents' (the event sequence number, matching `netname''s own chronological event order) - see {help nwdynam}.{p_end}"
			error 198
		}
		capture `__nwdynam_fpre' assert `__nwdynam_oppact' >= 1 & `__nwdynam_oppact' <= `nodes' if !missing(`__nwdynam_oppact')
		if _rc {
			di as error "{p 0 0 2}`__nwdynam_oppact' must be between 1 and `nodes' (the actor ID, matching `netname''s own actor order) - see {help nwdynam}.{p_end}"
			error 198
		}
		mata: `__nwdynam_oppmat' = J(`nevents', `nodes', 0)
		`__nwdynam_fpre' mata: st_view(__nwdynam_oppview=., ., "`__nwdynam_oppev' `__nwdynam_oppact'")
		mata: for (__nwdynam_oi=1; __nwdynam_oi<=rows(__nwdynam_oppview); __nwdynam_oi++) `__nwdynam_oppmat'[__nwdynam_oppview[__nwdynam_oi,1], __nwdynam_oppview[__nwdynam_oi,2]] = 1
		mata: mata drop __nwdynam_oppview __nwdynam_oi
		// as goldfish: the realised receiver of every event must be among its available actors (events in chronological
		// order; equal times keep the order of the event list)
		mata: __nwdynam_oord = order((`__nwdynam_events'[., 3], (1::rows(`__nwdynam_events'))), (1, 2))
		mata: __nwdynam_orcv = `__nwdynam_events'[__nwdynam_oord, 2]
		mata: __nwdynam_ook = J(rows(__nwdynam_orcv), 1, 0)
		mata: for (__nwdynam_oi = 1; __nwdynam_oi <= rows(__nwdynam_orcv); __nwdynam_oi++) __nwdynam_ook[__nwdynam_oi] = `__nwdynam_oppmat'[__nwdynam_oi, __nwdynam_orcv[__nwdynam_oi]]
		mata: st_local("__nwdynam_oppbad", strofreal(sum(__nwdynam_ook :== 0)))
		mata: st_local("__nwdynam_oppfirst", strofreal(min(select((1::rows(__nwdynam_ook)), __nwdynam_ook :== 0))))
		if `__nwdynam_oppbad' > 0 {
			mata: st_local("__nwdynam_opprcv", strofreal(__nwdynam_orcv[`__nwdynam_oppfirst']))
		}
		mata: mata drop __nwdynam_oord __nwdynam_orcv __nwdynam_ook __nwdynam_oi
		if `__nwdynam_oppbad' > 0 {
			_nwdynam_dropmata `__nwdynam_events' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2' `__nwdynam_oppmat' `__nwdynam_chg' `__nwdynam_comprows'
			di as error "{p 0 0 2}opportunities(): in `__nwdynam_oppbad' event(s) the realised receiver is not among the available actors (first: event `__nwdynam_oppfirst', receiver `__nwdynam_opprcv'); every event's receiver must be available.{p_end}"
			error 198
		}
	}
	else {
		mata: `__nwdynam_oppmat' = J(0, 0, 0)
	}

	// tie() ("cross-network effects, v1 scope", added 2026-09-02,
	// continuing the goldfish-comparison list after choice_coordination's
	// own effect list was completed) - goldfish's own `tie(network, ...)'
	// effect, `s(i,j,t,x) = I(x_ij>0)' (unweighted only in this v1),
	// identical in shape to `inertia()' but reading a SEPARATE,
	// EXOGENOUS network instead of the dependent one. v1 scope (a real,
	// disclosed limit, not silently narrow): `tie()' takes the NAME of
	// an already-declared, ORDINARY (non-event) network - resolved via
	// {help _nwsyntax}'s own `other()' prefix mechanism so its own
	// locals (`__nwdynam_tienet_netobj' etc.) never clash with the
	// primary network's own already-resolved locals. Genuinely
	// dynamically-evolving cross-network effects (an event-declared
	// SECOND network) are converted to a netterms() tie term earlier;
	// this block handles the STATIC network, verified
	// directly against real goldfish (see
	// unw_dynam.do's own header comment above the `tie_i' construction
	// for the exact formula and dev/dynam_unit18_tie_crosscheck.R/.do
	// for the real-goldfish verification). Valid for BOTH
	// submodel(choice) and submodel(choice_coordination) (added "batch
	// 16", 2026-09-02 - the coordination engine's own full-matrix
	// requirement made this the simplest possible addition: `S_tie' is
	// literally `tiemat' itself, no computation at all, verified
	// directly against real goldfish - see dev/dynam_unit19_
	// coordination_tie_crosscheck.R/.do), matching goldfish's own effect
	// table exactly; rejected under submodel(rate) (goldfish's own table
	// shows `tie' as choice/choice_coordination only).
	tempname __nwdynam_tiemat
	if "`tie'" != "" {
		if "`submodel'" == "rate" {
			di as error "{p 0 0 2}tie() only applies to submodel(choice) or submodel(choice_coordination) - see {help nwdynam}.{p_end}"
			error 198
		}
		_nwsyntax `tie', other(__nwdynam_tienet_) max(1)
		if "`__nwdynam_tienet_istemporal'" == "true" {
			di as error "{p 0 0 2}tie(`tie') is an event-type (eventtime()-declared) network - tie() only supports a STATIC exogenous network (a real, disclosed v1 scope limit; genuinely dynamically-evolving cross-network effects are not yet implemented) - see {help nwdynam}.{p_end}"
			error 198
		}
		if `__nwdynam_tienet_nodes' != `nodes' {
			di as error "{p 0 0 2}tie(`tie') has `__nwdynam_tienet_nodes' actors but `netname' has `nodes' - tie() requires the exogenous network to have exactly the same number of actors, in the same order (row/index correspondence, not label matching) - see {help nwdynam}.{p_end}"
			error 198
		}
		mata: `__nwdynam_tiemat' = (*(`__nwdynam_tienet_netobj'->get_matrix()) :> 0) :& (*(`__nwdynam_tienet_netobj'->get_matrix()) :!= .)
	}
	else {
		mata: `__nwdynam_tiemat' = J(0, 0, 0)
	}

	// Two-mode direction validation: goldfish's own architecture (and
	// this package's own contiguous mode-1-then-mode-2 actor-index
	// convention) assumes every event's sender is mode 1 (index <= n1)
	// and every receiver is mode 2 (index > n1) - never the reverse.
	// Checked directly against the real event data rather than assumed
	// true of any two-mode network a user might declare.
	if `__nwdynam_n1' > 0 {
		mata: st_numscalar("__nwdynam_baddir", sum((`__nwdynam_events'[.,1] :> `__nwdynam_n1') :| (`__nwdynam_events'[.,2] :<= `__nwdynam_n1')))
		local __nwdynam_baddir = __nwdynam_baddir
		scalar drop __nwdynam_baddir
		if `__nwdynam_baddir' > 0 {
			di as error "{p 0 0 2}nwdynam requires every event's sender to be mode 1 and receiver to be mode 2 for a two-mode (bipartite) network - `netname' has `__nwdynam_baddir' event(s) violating this (matching goldfish's own one-directional two-mode DyNAM architecture) - see {help nwdynam}.{p_end}"
			error 198
		}
	}

	// goldfish: an event of an actor who is absent at the time of the event is an error
	if `__nwdynam_ncomp' > 0 | "`present'" != "" {
		capture noisily mata: dynam_check_presence(`__nwdynam_comprows', `__nwdynam_events', `netobj'->get_nodenames(), `__nwdynam_chg')
		if _rc {
			local __nwdynam_crc = _rc
			_nwdynam_dropmata `__nwdynam_events' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2' `__nwdynam_chg' `__nwdynam_comprows'
			exit `__nwdynam_crc'
		}
	}
	mata: `__nwdynam_chg' = `__nwdynam_chg' \ `__nwdynam_comprows'
	capture mata: mata drop `__nwdynam_comprows'
	// scalar read when the event stream is prepared (option sequential); the fit runs inside capture so that a failure inside Mata
	// or the plugin still restores the user's frame (the native engine works in a frame of its own) and drops the helpers
	capture scalar drop __nwdynam_sequential
	scalar __nwdynam_sequential = ("`sequential'" != "")
	local __nwdynam_frame = c(frame)
	capture noisily {
		if "`submodel'" == "choice" {
			mata: DynamChoiceFitMulti(`__nwdynam_events', `nodes', (`__nwdynam_active'), `__nwdynam_same', `__nwdynam_diff', `__nwdynam_sim', `__nwdynam_alter', `__nwdynam_tertius', `__nwdynam_egoalterint1', `__nwdynam_egoalterint2', `__nwdynam_windowinertia', `__nwdynam_windowrecip', `__nwdynam_windowtertius', `__nwdynam_weightedinertia', `__nwdynam_weightedrecip', `__nwdynam_weightedindeg', `__nwdynam_weightedoutdeg', `__nwdynam_n1', `__nwdynam_oppmat', `__nwdynam_tiemat', `__nwdynam_chg', "__nwdynam_b", "__nwdynam_V", "__nwdynam_ll")
			local __nwdynam_title "Dynamic Network Actor Model, choice sub-model (conditional logit, MLE)"
			local __nwdynam_hdr "Dynamic Network Actor Model - choice sub-model (MLE)"
		}
		else if "`submodel'" == "rate" {
			if "`intercept'" != "" {
				mata: DynamRateInterceptFitMulti(`__nwdynam_events', `nodes', (`__nwdynam_active'), `__nwdynam_ego', `__nwdynam_tertius', `__nwdynam_windowindeg', `__nwdynam_windowoutdeg', `__nwdynam_windowtertius', `__nwdynam_weightedindeg', `__nwdynam_weightedoutdeg', `__nwdynam_n1', `__nwdynam_chg', "__nwdynam_b", "__nwdynam_V", "__nwdynam_ll")
				local __nwdynam_collabs "Intercept`__nwdynam_collabs'"
				local __nwdynam_title "Dynamic Network Actor Model, rate sub-model, WITH intercept (continuous-time competing-risks hazard, MLE)"
				local __nwdynam_hdr "Dynamic Network Actor Model - rate sub-model, with intercept (MLE)"
			}
			else {
				mata: DynamRateFitMulti(`__nwdynam_events', `nodes', (`__nwdynam_active'), `__nwdynam_ego', `__nwdynam_tertius', `__nwdynam_windowindeg', `__nwdynam_windowoutdeg', `__nwdynam_windowtertius', `__nwdynam_weightedindeg', `__nwdynam_weightedoutdeg', `__nwdynam_n1', `__nwdynam_chg', "__nwdynam_b", "__nwdynam_V", "__nwdynam_ll")
				local __nwdynam_title "Dynamic Network Actor Model, rate sub-model, no intercept (ordinal partial likelihood, MLE)"
				local __nwdynam_hdr "Dynamic Network Actor Model - rate sub-model, no intercept (MLE)"
			}
		}
		else {
			mata: DynamCoordFitMulti(`__nwdynam_events', `nodes', (`__nwdynam_active'), `__nwdynam_same', `__nwdynam_diff', `__nwdynam_sim', `__nwdynam_alter', `__nwdynam_tertius', `__nwdynam_egoalterint1', `__nwdynam_egoalterint2', `__nwdynam_windowinertia', `__nwdynam_windowindeg', `__nwdynam_windowtertius', `__nwdynam_weightedinertia', `__nwdynam_weightedindeg', `__nwdynam_tiemat', `__nwdynam_chg', "__nwdynam_b", "__nwdynam_V", "__nwdynam_ll")
			local __nwdynam_title "Dynamic Network Actor Model, choice_coordination sub-model (multinomial-multinomial, MLE)"
			local __nwdynam_hdr "Dynamic Network Actor Model - choice_coordination sub-model (MLE)"
		}
	}
	local __nwdynam_rc = _rc
	capture scalar drop __nwdynam_sequential
	capture mata: mata drop __nwdynam_xblob
	if `__nwdynam_rc' {
		if "`c(frame)'" != "`__nwdynam_frame'" capture frame change `__nwdynam_frame'
		capture frame drop __nwdynam_native
		_nwdynam_dropmata `__nwdynam_events' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2' `__nwdynam_oppmat' `__nwdynam_tiemat' `__nwdynam_chg' `__nwdynam_comprows'
		capture matrix drop __nwdynam_b
		capture matrix drop __nwdynam_V
		capture scalar drop __nwdynam_ll
		exit `__nwdynam_rc'
	}
	// the fit's data stays in frame __nwdynam_fit for estat; add the actor labels and the event order to it (nothing of this is part of the estimation)
	capture noisily mata: st_local("__nwdynam_fitid", dynam_fit_meta(`netobj'->get_nodenames(), `__nwdynam_events'))
	if _rc local __nwdynam_fitid ""
	mata: mata drop `__nwdynam_events' `__nwdynam_same' `__nwdynam_diff' `__nwdynam_sim' `__nwdynam_ego' `__nwdynam_alter' `__nwdynam_tertius' `__nwdynam_egoalterint1' `__nwdynam_egoalterint2' `__nwdynam_oppmat' `__nwdynam_tiemat' `__nwdynam_chg'

	tempname b V
	matrix `b' = __nwdynam_b
	matrix `V' = __nwdynam_V
	matrix colnames `b' = `__nwdynam_collabs'
	matrix colnames `V' = `__nwdynam_collabs'
	matrix rownames `V' = `__nwdynam_collabs'
	local __nwdynam_ll = __nwdynam_ll
	matrix drop __nwdynam_b __nwdynam_V
	scalar drop __nwdynam_ll

	ereturn post `b' `V', obs(`nevents')
	ereturn local cmd "nwdynam"
	ereturn local version "1.9.1"
	ereturn local title "`__nwdynam_title'"
	ereturn local depvar "`netname'"
	ereturn local submodel "`submodel'"
	ereturn local effects "`=trim("`__nwdynam_collabs'")'"
	ereturn scalar N = `nevents'
	ereturn scalar nodes = `nodes'
	ereturn scalar ll = `__nwdynam_ll'
	ereturn scalar ntimes = `__nwdynam_ntimes'
	if "`present'`compchange'" != "" {
		ereturn scalar ncompchange = `__nwdynam_ncomp'
		ereturn local present "`present'"
		ereturn local compchange "`compchange'"
	}
	ereturn local fitid "`__nwdynam_fitid'"
	ereturn local estat_cmd "nwdynam_estat"
	ereturn local ties = cond(`__nwdynam_ntimes' < `nevents', cond("`sequential'" != "", "sequential", "simultaneous"), "none")

	ereturn local header "`__nwdynam_hdr'"
	// effects that are not identified: no standard error, a zero variance (flat likelihood) or an enormous standard
	// error (separation, collinear effects)
	local __nwdynam_bad ""
	local __nwdynam_k = colsof(e(b))
	forvalues __nwdynam_j = 1/`__nwdynam_k' {
		local __nwdynam_v = el(e(V), `__nwdynam_j', `__nwdynam_j')
		if missing(`__nwdynam_v') | `__nwdynam_v' <= 0 {
			local __nwdynam_bad "`__nwdynam_bad' `: word `__nwdynam_j' of `__nwdynam_collabs''"
		}
		else if sqrt(`__nwdynam_v') > 100 {
			local __nwdynam_bad "`__nwdynam_bad' `: word `__nwdynam_j' of `__nwdynam_collabs''"
		}
	}
	ereturn local unidentified "`=trim("`__nwdynam_bad'")'"
	_nwdynam_display, `__tabopts'
end

capture program drop _nwdynam_display
program _nwdynam_display
	syntax [, *]
	di as text "{hline 60}"
	di as text "`e(header)'"
	di as text "Network: " as result "`e(depvar)'" _col(40) as text "Actors: " as result e(nodes)
	di as text "Events: " as result e(N) _col(40) as text "Log likelihood: " as result %9.4f e(ll)
	if e(ntimes) < e(N) {
		if "`e(ties)'" == "sequential" di as text "{p 0 0 2}Note: " as result e(N) - e(ntimes) as text " events share a time point; processed one after another (option sequential).{p_end}"
		else di as text "{p 0 0 2}Note: " as result e(N) - e(ntimes) as text " events share a time point with an earlier event; events with equal times are simultaneous (as in goldfish).{p_end}"
	}
	if "`e(present)'`e(compchange)'" != "" {
		di as text "{p 0 0 2}Composition change: " as result e(ncompchange) as text " join/leave events" ///
			cond("`e(compchange)'" != "", " (" + "`e(compchange)'" + ")", "") ///
			cond("`e(present)'" != "", "; presence at the start: " + "`e(present)'", "") as text ".{p_end}"
	}
	di as text "{hline 60}"
	_nwdynam_tabopts __t, `options'
	ereturn display, `__t'
	if "`e(unidentified)'" != "" {
		di as text "{p 0 6 2}Note: no usable standard error for " as result "`e(unidentified)'" as text ".{p_end}"
		di as text "{p 6 6 2}These effects are probably not identified (too few events, a perfectly predicted pattern, or effects"   ///
			" that are linear functions of each other). Do not interpret their estimates; drop or recode the effect.{p_end}"
	}
end

capture program drop _nwdynam_dropmata
program _nwdynam_dropmata
	foreach __m of local 0 {
		capture mata: mata drop `__m'
	}
end

// the options of ereturn display for the coefficient table. The coefficient names are shown in full: the
// confidence intervals are left out (noci) at a linesize below 100, where names of up to 32 characters do not fit next
// to them. ci, or level(#) given explicitly, shows the intervals (long names are then cut, e.g. nodematch_o~e); noci,
// nopvalues, cformat() and the other display options of Stata's estimation commands are passed on. With check, only
// the options are checked (before a fit).
// Use: _nwdynam_tabopts <local> [check] [, options]; the options for ereturn display are returned in <local>.
capture program drop _nwdynam_tabopts
program _nwdynam_tabopts
	version 14
	gettoken __mac 0 : 0, parse(" ,")
	gettoken __chk : 0, parse(" ,")
	if `"`__chk'"' == "check" gettoken __chk 0 : 0, parse(" ,")
	else local __chk
	syntax [, Level(string) *]
	local __lev `"`level'"'
	// ci is picked out by hand: declared to syntax, ci would also take noci (as its negation)
	local ci
	local __rest
	while `"`options'"' != "" {
		gettoken __o options : options, bind
		if strlower(`"`__o'"') == "ci" local ci ci
		else local __rest `"`__rest' `__o'"'
	}
	_get_diopts __di, `__rest'
	local __noci : list posof "noci" in __di
	local __nopv : list posof "nopvalues" in __di
	if "`ci'" != "" & `__noci' {
		di as err "only one of options ci or noci is allowed"
		exit 198
	}
	local __show = ("`ci'" != "" | `"`__lev'"' != "")
	if `"`__lev'"' != "" {
		local 0 , level(`__lev')
		syntax [, Level(cilevel)]
		local __lev `level'
	}
	else local __lev = c(level)
	if "`__chk'" != "" exit
	local __add
	if !`__show' & !`__noci' & !`__nopv' {
		// full names: at a linesize below 100 (or with lstretch off) the names of up to 32 characters do not fit next
		// to the intervals; from 100 on both fit
		if c(linesize) < 100 | c(lstretch) == "off" | `: list posof "nolstretch" in __di' local __add noci
	}
	c_local `__mac' `"level(`__lev') `__di' `__add'"'
end
