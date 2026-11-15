/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
 * nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	dynam_sim.c -- native (C) DyNAM likelihood/gradient kernel for
	nwcommands (nwdynam).

	See the development notes (scope/status) and unw_dynam.do's own
	header comment (Mata reference implementation, built and
	CORRECTNESS-VERIFIED against a real reference R implementation
	FIRST - dev/dynam_unit1_choice_crosscheck.do,
	dev/dynam_unit2_rate_crosscheck.do - matching this project's own
	established sequencing discipline, ergm_mcmc.c's own header: get the
	Mata reference right and crosschecked before writing
	performance-critical C, since a fast wrong answer is worse than a
	slow right one).

	WHAT THIS ACCELERATES: unlike ergm_mcmc.c (a full from-scratch MCMC
	sampler, since Mata's own optimize() has no equivalent) or
	saom_sim.c (a full ministep simulation loop), nwdynam fits via plain
	MLE - Mata's own optimize() (BFGS, with a Nelder-Mead fallback)
	already handles the outer optimization robustly and is NOT
	reimplemented here. This plugin accelerates only the HOT INNER LOOP:
	one log-likelihood + analytic-gradient evaluation at a given theta,
	called repeatedly (once per outer optimizer iteration - typically a
	few dozen calls per fit, not thousands) from Mata's own evaluator
	callback. Each such call still does the SAME O(n) per-event work the
	Mata engine does (dynam_choice_loglik_grad_unit1()/
	dynam_rate_loglik_grad_unit1() in unw_dynam.do) - the win is
	removing Mata's own per-operation interpreter overhead on that
	tight, compounding loop, the identical rationale nwgraph.c's own
	header already gives for betweenness centrality.

	ONE REAL ALGORITHMIC IMPROVEMENT OVER THE MATA REFERENCE, NOT JUST A
	COMPILED PORT: unw_dynam.do's own dynam_choice_loglik_grad_unit1()
	recomputes `indeg_i = colsum(tie)` - an O(n) reduction over an n x n
	MATRIX, i.e. genuinely O(n^2) work - AT EVERY EVENT, for both the
	choice sub-model's own indeg-alter effect and (implicitly, via the
	analogous call) the rate sub-model's indeg/outdeg-ego effects. Since
	`tie` only ever grows (a 0 cell flips to 1 exactly once, on the
	first tie between that ordered pair; a repeat contact leaves it at
	1), each actor's own in-degree and out-degree can instead be
	maintained as a genuinely O(1)-per-update running count (`indeg[]`/
	`outdeg[]` below), collapsing an O(n^2 * nevents) total cost down to
	O(n * nevents) - the same asymptotic class the row/column inertia/
	recip lookups already need per event regardless. This is a real,
	disclosed finding from porting the reference algorithm to C, not
	silently different behavior: both implementations compute
	MATHEMATICALLY IDENTICAL quantities (indeg[j] here always equals
	colsum(tie)[j] there, by construction - tie is monotonically
	non-decreasing and indeg[j] is incremented exactly once per 0->1
	flip in column j), so the CERTIFICATION CONTRACT below still holds
	exactly, not approximately.

	CERTIFICATION CONTRACT: like nwgraph.c's betweenness kernel (and
	unlike ergm_mcmc.c's stochastic MCMC sampler), this is a
	DETERMINISTIC, exact numerical computation - the native and Mata
	implementations of the same log-likelihood/gradient formula on the
	same event data and the same theta must agree EXACTLY (up to
	ordinary floating-point summation-order effects from the
	incremental-vs-colsum reformulation above, not a statistical
	tolerance) - see cscripts/test_nwdynam_native.do.

	WIRE PROTOCOL: one `plugin call' per OUTER OPTIMIZER ITERATION (not
	per event) - see unw_dynam.do's own DynamChoiceEvalNative()/
	DynamRateEvalNative(). Two input frame variables, set up ONCE before
	the optimizer loop begins and read fresh on every call (a plugin
	invocation carries no state across calls, same as nwgraph.c's own
	convention - the per-call CSR/tie-matrix rebuild is a real, accepted
	cost, small next to the loop it feeds and tiny next to Mata's own
	per-operation interpreter overhead this whole native path exists to
	avoid): v1=sender, v2=receiver (nevents rows, 1-indexed actor IDs,
	in chronological order - the CALLER is responsible for sorting, this
	file does not re-sort). One output frame variable (v3), used to
	return a small fixed-shape result: row 1 = log-likelihood, rows
	2..(nparams+1) = the analytic gradient, in the same effect order
	unw_dynam.do's own theta vectors use (choice: inertia, recip, indeg;
	rate: indeg, outdeg). Args string: "algcode n nevents theta1
	[theta2 [theta3]]" (theta values passed as text, since they are few
	in number and change on every call - cheaper than a frame column
	round-trip for 2-3 numbers).
*/

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Plugin codes 11 and 12 = codes 1 and 2 with SIMULTANEOUS events (goldfish): a fourth variable v4 is 1 at the last
   event of each group of equal times, and the tie state is updated only after the last event of a group.
   Adding 100 to a code (101, 102, 111, 112) passes the VALUE of each event in a fifth variable v5 (v4 is then always
   passed): as in goldfish, the value is the increment of the tie, a tie exists while the summed value is above 0, and
   indeg/outdeg add the sign of each tie value (-1 below 0).
   Without it every event has the value 1 (a tie exists after the first event of a pair). */
#define ALG_DYNAM_CHOICE 1
#define ALG_DYNAM_RATE 2
#define ALG_DYNAM_RATEINT 3

/* ===================================================================
   Choice sub-model: inertia + recip + indeg(alter). Direct port of
   dynam_choice_loglik_grad_unit1() (unw_dynam.do), with indeg tracked
   incrementally (see this file's own header) instead of a fresh
   colsum() every event.
   =================================================================== */

static int dynam_choice_loglik_grad(long n, long nevents, long *sender, long *receiver, char *gend, double *val,
		double theta1, double theta2, double theta3,
		double *out_ll, double out_grad[3]) {
	char *tie = NULL;      /* n x n, 0/1, row-major: tie[(s-1)*n + (j-1)] */
	double *xval = NULL;   /* n x n summed event values (only with val) */
	double *indeg = NULL;  /* n, running in-degree (distinct senders so far) */
	double *lp = NULL, *P = NULL;
	long i, j, k, s, r, gs = 0;
	double ll = 0.0, grad1 = 0.0, grad2 = 0.0, grad3 = 0.0;
	double mx, sumP, lrsum, inertia_r, recip_r, indeg_r;
	int rc = 0;

	tie = (char *)calloc((size_t)n * (size_t)n, sizeof(char));
	if (val) xval = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
	indeg = (double *)calloc((size_t)n, sizeof(double));
	lp = (double *)malloc((size_t)n * sizeof(double));
	P = (double *)malloc((size_t)n * sizeof(double));
	if (!tie || (val && !xval) || !indeg || !lp || !P) { rc = 909; goto cleanup; }

	for (i = 0; i < nevents; i++) {
		s = sender[i];
		r = receiver[i];

		mx = -1e300;
		for (j = 1; j <= n; j++) {
			if (j == s) { lp[j - 1] = -1e300; continue; }
			{
				double inertia_j = tie[(s - 1) * n + (j - 1)] ? 1.0 : 0.0;
				double recip_j = tie[(j - 1) * n + (s - 1)] ? 1.0 : 0.0;
				double indeg_j = indeg[j - 1];
				lp[j - 1] = theta1 * inertia_j + theta2 * recip_j + theta3 * indeg_j;
			}
			if (lp[j - 1] > mx) mx = lp[j - 1];
		}

		sumP = 0.0;
		for (j = 1; j <= n; j++) {
			P[j - 1] = exp(lp[j - 1] - mx);
			sumP += P[j - 1];
		}
		lrsum = mx + log(sumP);
		ll += (lp[r - 1] - lrsum);

		{
			double sum_inertia_P = 0.0, sum_recip_P = 0.0, sum_indeg_P = 0.0;
			for (j = 1; j <= n; j++) {
				double pj = P[j - 1] / sumP;
				double inertia_j = (j == s) ? 0.0 : (tie[(s - 1) * n + (j - 1)] ? 1.0 : 0.0);
				double recip_j = (j == s) ? 0.0 : (tie[(j - 1) * n + (s - 1)] ? 1.0 : 0.0);
				double indeg_j = (j == s) ? 0.0 : indeg[j - 1];
				sum_inertia_P += inertia_j * pj;
				sum_recip_P += recip_j * pj;
				sum_indeg_P += indeg_j * pj;
			}
			inertia_r = tie[(s - 1) * n + (r - 1)] ? 1.0 : 0.0;
			recip_r = tie[(r - 1) * n + (s - 1)] ? 1.0 : 0.0;
			indeg_r = indeg[r - 1];
			grad1 += (inertia_r - sum_inertia_P);
			grad2 += (recip_r - sum_recip_P);
			grad3 += (indeg_r - sum_indeg_P);
		}

		/* simultaneous events (same time): the state is updated only after the last event of the group */
		if (!gend || gend[i]) {
			for (k = gs; k <= i; k++) {
				long sk = sender[k], rk = receiver[k];
				if (val) {
					long c = (sk - 1) * n + (rk - 1);
					char now;
					/* indeg adds the sign of each tie value (goldfish: -1 below 0), inertia/recip test > 0 */
					double sold = (xval[c] > 0.0) - (xval[c] < 0.0), snew;
					xval[c] += val[k];
					now = (xval[c] > 0.0);
					snew = (xval[c] > 0.0) - (xval[c] < 0.0);
					indeg[rk - 1] += snew - sold;
					tie[c] = now;
				}
				else if (!tie[(sk - 1) * n + (rk - 1)]) {
					tie[(sk - 1) * n + (rk - 1)] = 1;
					indeg[rk - 1] += 1.0;
				}
			}
			gs = i + 1;
		}
	}

	*out_ll = ll;
	out_grad[0] = grad1;
	out_grad[1] = grad2;
	out_grad[2] = grad3;

cleanup:
	free(tie);
	free(xval);
	free(indeg);
	free(lp);
	free(P);
	return rc;
}

/* ===================================================================
   Rate sub-model (no intercept): indeg + outdeg, both ego type. Direct
   port of dynam_rate_loglik_grad_unit1() (unw_dynam.do) - risk set is
   ALL n actors (no self-exclusion; the realized sender IS one of the n
   candidates, unlike the choice sub-model's n-1). Both indeg and
   outdeg tracked incrementally throughout, same reasoning as the
   choice kernel above.
   =================================================================== */

static int dynam_rate_loglik_grad(long n, long nevents, long *sender, long *receiver, char *gend, double *val,
		double theta1, double theta2,
		double *out_ll, double out_grad[2]) {
	char *tie = NULL;
	double *xval = NULL;
	double *indeg = NULL, *outdeg = NULL;
	double *lp = NULL, *P = NULL;
	long i, j, k, s, r, gs = 0;
	double ll = 0.0, grad1 = 0.0, grad2 = 0.0;
	double mx, sumP, lrsum;
	int rc = 0;

	tie = (char *)calloc((size_t)n * (size_t)n, sizeof(char));
	if (val) xval = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
	indeg = (double *)calloc((size_t)n, sizeof(double));
	outdeg = (double *)calloc((size_t)n, sizeof(double));
	lp = (double *)malloc((size_t)n * sizeof(double));
	P = (double *)malloc((size_t)n * sizeof(double));
	if (!tie || (val && !xval) || !indeg || !outdeg || !lp || !P) { rc = 909; goto cleanup; }

	for (i = 0; i < nevents; i++) {
		s = sender[i];
		r = receiver[i];

		mx = -1e300;
		for (j = 1; j <= n; j++) {
			lp[j - 1] = theta1 * indeg[j - 1] + theta2 * outdeg[j - 1];
			if (lp[j - 1] > mx) mx = lp[j - 1];
		}

		sumP = 0.0;
		for (j = 1; j <= n; j++) {
			P[j - 1] = exp(lp[j - 1] - mx);
			sumP += P[j - 1];
		}
		lrsum = mx + log(sumP);
		ll += (lp[s - 1] - lrsum);

		{
			double sum_indeg_P = 0.0, sum_outdeg_P = 0.0;
			for (j = 1; j <= n; j++) {
				double pj = P[j - 1] / sumP;
				sum_indeg_P += indeg[j - 1] * pj;
				sum_outdeg_P += outdeg[j - 1] * pj;
			}
			grad1 += (indeg[s - 1] - sum_indeg_P);
			grad2 += (outdeg[s - 1] - sum_outdeg_P);
		}

		if (!gend || gend[i]) {
			for (k = gs; k <= i; k++) {
				long sk = sender[k], rk = receiver[k];
				if (val) {
					long c = (sk - 1) * n + (rk - 1);
					char now;
					/* indeg/outdeg add the sign of each tie value (goldfish: -1 below 0) */
					double sold = (xval[c] > 0.0) - (xval[c] < 0.0), snew;
					xval[c] += val[k];
					now = (xval[c] > 0.0);
					snew = (xval[c] > 0.0) - (xval[c] < 0.0);
					indeg[rk - 1] += snew - sold;
					outdeg[sk - 1] += snew - sold;
					tie[c] = now;
				}
				else if (!tie[(sk - 1) * n + (rk - 1)]) {
					tie[(sk - 1) * n + (rk - 1)] = 1;
					indeg[rk - 1] += 1.0;
					outdeg[sk - 1] += 1.0;
				}
			}
			gs = i + 1;
		}
	}

	*out_ll = ll;
	out_grad[0] = grad1;
	out_grad[1] = grad2;

cleanup:
	free(tie);
	free(xval);
	free(indeg);
	free(outdeg);
	free(lp);
	free(P);
	return rc;
}


/* ===================================================================
   Rate sub-model WITH intercept (real-time hazard), effects indeg and outdeg (nwdynam 1.4.0). Port of
   dynam_rateint_loglik_grad_multi() (unw_dynam.do) for the unweighted, unwindowed, one-mode case. Log likelihood
     sum_i [ lp(s_i) - w_i * sum_j exp(lp_j) ],  lp_j = th0 + th1 * indeg_j + th2 * outdeg_j,
   w_i = time since the previous event (0 for the first event and inside a group of equal times). Degrees are kept
   up to date in O(1) per tie update, so one evaluation costs O(n * nevents) time and O(n^2) bytes, with no cache of
   per-event statistics. Returns the log likelihood, the gradient (3) and the full Hessian (3 x 3, upper triangle
   by rows: 00 01 02 11 12 22). If lp overflows exp() the log likelihood is returned missing, so that the optimizer
   shortens its step.
   =================================================================== */

static int dynam_rateint_loglik_grad(long n, long nevents, long *sender, long *receiver, char *gend, double *val,
		double *w, double th[3], int *overflow, double *out_ll, double out_grad[3], double out_h[6]) {
	char *tie = NULL;
	double *xval = NULL;
	double *indeg = NULL, *outdeg = NULL, *e = NULL;
	long i, j, k, s, gs = 0;
	double ll = 0.0, g[3] = {0, 0, 0}, h[6] = {0, 0, 0, 0, 0, 0};
	int rc = 0;

	*overflow = 0;
	tie = (char *)calloc((size_t)n * (size_t)n, sizeof(char));
	if (val) xval = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
	indeg = (double *)calloc((size_t)n, sizeof(double));
	outdeg = (double *)calloc((size_t)n, sizeof(double));
	e = (double *)malloc((size_t)n * sizeof(double));
	if (!tie || (val && !xval) || !indeg || !outdeg || !e) { rc = 909; goto cleanup; }

	for (i = 0; i < nevents; i++) {
		double wi = w[i], lps, S = 0.0, Sx1 = 0.0, Sx2 = 0.0, Sx11 = 0.0, Sx12 = 0.0, Sx22 = 0.0;
		s = sender[i] - 1;
		lps = th[0] + th[1] * indeg[s] + th[2] * outdeg[s];
		if (lps > 700.0) { *overflow = 1; goto cleanup; }
		ll += lps;
		g[0] += 1.0;
		g[1] += indeg[s];
		g[2] += outdeg[s];
		if (wi != 0.0) {
			for (j = 0; j < n; j++) {
				double lp = th[0] + th[1] * indeg[j] + th[2] * outdeg[j];
				double x1 = indeg[j], x2 = outdeg[j], ej;
				if (lp > 700.0) { *overflow = 1; goto cleanup; }
				ej = exp(lp);
				S += ej;
				Sx1 += ej * x1;
				Sx2 += ej * x2;
				Sx11 += ej * x1 * x1;
				Sx12 += ej * x1 * x2;
				Sx22 += ej * x2 * x2;
			}
			ll -= wi * S;
			g[0] -= wi * S;
			g[1] -= wi * Sx1;
			g[2] -= wi * Sx2;
			h[0] -= wi * S;
			h[1] -= wi * Sx1;
			h[2] -= wi * Sx2;
			h[3] -= wi * Sx11;
			h[4] -= wi * Sx12;
			h[5] -= wi * Sx22;
		}

		if (!gend || gend[i]) {
			for (k = gs; k <= i; k++) {
				long sk = sender[k], rk = receiver[k];
				if (val) {
					long c = (sk - 1) * n + (rk - 1);
					double sold = (xval[c] > 0.0) - (xval[c] < 0.0), snew;
					xval[c] += val[k];
					snew = (xval[c] > 0.0) - (xval[c] < 0.0);
					indeg[rk - 1] += snew - sold;
					outdeg[sk - 1] += snew - sold;
					tie[c] = (xval[c] > 0.0);
				}
				else if (!tie[(sk - 1) * n + (rk - 1)]) {
					tie[(sk - 1) * n + (rk - 1)] = 1;
					indeg[rk - 1] += 1.0;
					outdeg[sk - 1] += 1.0;
				}
			}
			gs = i + 1;
		}
	}
	*out_ll = ll;
	for (k = 0; k < 3; k++) out_grad[k] = g[k];
	for (k = 0; k < 6; k++) out_h[k] = h[k];

cleanup:
	free(tie);
	free(xval);
	free(indeg);
	free(outdeg);
	free(e);
	return rc;
}

/* ===================================================================
   Plugin entry point
   =================================================================== */

static double next_double(void) {
	return atof(strtok(NULL, " \t"));
}

STDLL stata_call(int argc, char *argv[]) {
	char *argbuf;
	long algcode, n, nevents, i;
	int grouped = 0, valued = 0;
	char *gend = NULL;
	double *val = NULL;
	double theta1 = 0.0, theta2 = 0.0, theta3 = 0.0;
	long *sender = NULL, *receiver = NULL;
	double ll = 0.0, grad[3], hess[6], th3[3];
	double *wgap = NULL;
	int rc = 0, overflow = 0;

	if (argc < 1) { SF_error("dynam_sim: missing argument string\n"); return(198); }

	argbuf = (char *)malloc(strlen(argv[0]) + 1);
	strcpy(argbuf, argv[0]);
	algcode = (long)atof(strtok(argbuf, " \t"));
	/* codes 11/12: the same kernels with simultaneous events; v4 = 1 at the last event of each group of equal times */
	if (algcode > 100) { valued = 1; algcode -= 100; }
	if (algcode > 10) { grouped = 1; algcode -= 10; }
	n       = (long)next_double();
	nevents = (long)next_double();
	theta1  = next_double();
	if (algcode == ALG_DYNAM_CHOICE) {
		theta2 = next_double();
		theta3 = next_double();
	}
	else if (algcode == ALG_DYNAM_RATE) {
		theta2 = next_double();
	}
	else if (algcode == ALG_DYNAM_RATEINT) {
		theta2 = next_double();
		theta3 = next_double();
	}
	free(argbuf);

	if (nevents > 0) {
		sender = (long *)malloc((size_t)nevents * sizeof(long));
		receiver = (long *)malloc((size_t)nevents * sizeof(long));
		if (grouped) gend = (char *)malloc((size_t)nevents);
		if (valued) val = (double *)malloc((size_t)nevents * sizeof(double));
		if (algcode == ALG_DYNAM_RATEINT) wgap = (double *)malloc((size_t)nevents * sizeof(double));
		for (i = 0; i < nevents; i++) {
			ST_double v;
			SF_vdata(1, i + 1, &v); sender[i] = (long)v;
			SF_vdata(2, i + 1, &v); receiver[i] = (long)v;
			if (grouped) { SF_vdata(4, i + 1, &v); gend[i] = (v != 0.0); }
			if (valued) { SF_vdata(5, i + 1, &v); val[i] = v; }
			if (wgap) { SF_vdata(6, i + 1, &v); wgap[i] = v; }
		}
	}

	switch (algcode) {
		case ALG_DYNAM_CHOICE:
			rc = dynam_choice_loglik_grad(n, nevents, sender, receiver, gend, val, theta1, theta2, theta3, &ll, grad);
			if (rc == 0) {
				SF_vstore(3, 1, ll);
				SF_vstore(3, 2, grad[0]);
				SF_vstore(3, 3, grad[1]);
				SF_vstore(3, 4, grad[2]);
			}
			break;
		case ALG_DYNAM_RATE:
			rc = dynam_rate_loglik_grad(n, nevents, sender, receiver, gend, val, theta1, theta2, &ll, grad);
			if (rc == 0) {
				SF_vstore(3, 1, ll);
				SF_vstore(3, 2, grad[0]);
				SF_vstore(3, 3, grad[1]);
			}
			break;
		case ALG_DYNAM_RATEINT:
			th3[0] = theta1; th3[1] = theta2; th3[2] = theta3;
			rc = dynam_rateint_loglik_grad(n, nevents, sender, receiver, gend, val, wgap, th3, &overflow, &ll, grad, hess);
			if (rc == 0) {
				int q;
				SF_vstore(3, 1, overflow ? SV_missval : ll);
				for (q = 0; q < 3; q++) SF_vstore(3, 2 + q, overflow ? SV_missval : grad[q]);
				for (q = 0; q < 6; q++) SF_vstore(3, 5 + q, overflow ? SV_missval : hess[q]);
			}
			break;
		default:
			SF_error("dynam_sim: unknown algorithm code\n");
			rc = 198;
	}

	free(sender);
	free(receiver);
	free(gend);
	free(val);
	free(wgap);
	return rc;
}
