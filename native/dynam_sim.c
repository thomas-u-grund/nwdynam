/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
 * nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	dynam_sim.c -- native (C) log likelihood, gradient and Hessian of every nwdynam model
	(choice, rate, rate with intercept, choice_coordination).

	Every fit of nwdynam runs here: Stata's Mata optimize() drives the outer loop and calls this plugin once per
	evaluation at a given theta. The plugin carries no state between calls. The model (events, covariates, windows,
	masks) is written ONCE per fit as a block of numbers into the first variable of a frame (the "blob", layout below);
	each call passes only the length of the blob, the parameter count and theta in the argument string and gets the
	result in the second variable: row 1 = log likelihood (missing when exp() overflows), rows 2..P+1 = gradient, then
	the P x P Hessian by rows (only when todo = 2; choice_coordination has no Hessian).

	Argument string: "<code> <L> <todo> <P> <theta_1> ... <theta_P>"
	  code 20 = evaluate; code 99 = report the protocol number (row 1 of variable 2). Version 1 plugins (nwdynam 1.4.0
	  and earlier) know the codes 1, 2, 3 only and fail on 99, which is how the Mata glue recognises a plugin that is too old.

	Codes 30 and 40 do not change the blob or the codes above.
	  code 30 = per-event results after a fit (goldfish's returnIntervalLogL and returnEventProbabilities):
	  "30 <L> <flags> <P> <theta...> [<nrec>]". flags bit 1 = write the rows (else only count the records; with bit 1 the last token is
	  the record count of the counting call), bit 2 = write the probabilities too. Output: variable 2 rows 1..4 = log likelihood,
	  number of records, number of dependent events, values per probability block; the rows are written to variable 1 BEHIND the blob
	  (row L + 1 on): the log likelihood of every record (nrec rows), then the event number of every record (1-based position in the
	  time-sorted list; 0 = a right-censored interval of the rate model with intercept, which goldfish also lists), then the time at
	  the end of the record, then the probabilities (blocks of n values per event; for choice_coordination n(n-1)/2 values of the
	  pairs a < b; a missing value = not in the risk set; rate model with intercept: the linear predictor minus the exposure,
	  goldfish's value there).
	  code 40 = series tools on the per-event log likelihood (goldfish's examineOutliers and examineChangepoints, the latter after
	  the R package 'changepoint'): "40 <L> <nd> <kind> <method> <moment> <minseglen> <parameter> <window>"; variable 1 behind the
	  blob: the series (nd rows), the times (nd rows), then the result indices; variable 2 row 1 = number of indices, row 2 = error code.

	Blob (doubles; L = total length):
	  [0] protocol  [1] model (1 choice, 2 rate, 3 rate with intercept, 4 coordination)  [2] n  [3] N  [4] n1 (mode-1
	  actors of a two-mode network, 0 = one-mode)  [5] M (number of effect slots)  [6] opportunities given  [7] tie() given
	  [8..11]  windows inertia, recip, indeg, outdeg (>= 1e300 = none)
	  [12..15] weighted flags inertia, recip, indeg, outdeg
	  [16..35] active flags of the effect slots (choice 17, rate 5, coordination 12)
	  [36] C (number of change events of covariates)  [37] Cp (number of composition changes)  [38] initial presence given
	  [39] window of tertius (>= 1e300 = none)  [40] length of the block of exogenous networks and extra terms (0 = none)
	  then events: N rows of (sender, receiver, time, value, last-of-group flag), sorted by time (the caller sorts)
	  of n values (same, diff, sim, alter/ego, tertius, egoalterint1, egoalterint2): the values at the start of the log,
	  then the N x n opportunity matrix (only when [6] = 1), then the n x n matrix of tie() (only when [7] = 1).

	Composition change (as goldfish): the actors have a presence flag; an absent actor is not in the risk set (rate: the senders
	of the hazard and of the normaliser; choice: the candidate receivers; coordination: both members of a pair), its statistics
	are still computed. The blob carries, after the change events, Cp composition changes (actor, time, new presence 0/1) sorted
	by time and, when [38] = 1, the n values of the presence at the start (else everybody is present). A composition change takes
	effect for every event row whose time is LARGER OR EQUAL to its time (goldfish, unlike covariate changes, applies the changes
	of a time point BEFORE the events of that time point); in the rate model with intercept the whole period that ends at an event
	(or at a covariate change) uses the presence valid at the end of the period. The caller guarantees that the sender and
	receiver of every event are present; the plugin refuses an event otherwise.

	Time-varying covariates (as goldfish): a change event sets the value of a covariate row for one actor; it takes effect for
	the events whose time is LARGER than its time (an event at exactly the same time is evaluated with the old value: in goldfish
	the dependent events of a time point are processed before the other event lists). Changes before the first event and after
	the last one are allowed. The caller turns increments into new values. In the rate model with intercept the time between two
	events is split at every change, because the hazard of the intervals depends on the covariates of that moment.

	The statistics are those of the Mata reference evaluators (cscripts/ref_dynam_mata.do), which agree with goldfish.
	All counts are kept up to date in O(1) or O(n) per tie change instead of being recomputed at every event.
*/

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

typedef int64_t nwint;

#define DYNAM_PROTOCOL 5
#define NOWIN 1e300
#define HDR 41
#define XMAXL 16
#define XMAXT 32
#define XMAXA 16
#define XT_TIE 1
#define XT_TRANS 2
#define XT_CYCLE 3
#define XT_CSEND 4
#define XT_CRECV 5
#define XT_MTRANS 6
#define XT_MCYCLE 7
#define XT_MCSEND 8
#define XT_MCRECV 9
#define XT_TDIFF 10
#define XT_FOUR 11
#define M_CHOICE 1
#define M_RATE 2
#define M_RATEINT 3
#define M_COORD 4

typedef struct { int hasmat, sym; nwint ne; const double *ev; const double *mat; } XLayer;
typedef struct { int type, la, lb, ign, attr; double win, wgt; } XTermSpec;

typedef struct {
	int model;
	nwint n, N, n1;
	int M, hasopp, hastie;
	double win[5], wgt[4];      /* windows: inertia, recip, indeg, outdeg, tertius */
	double act[20];
	double *ev;        /* N x 5 */
	double *cov[7];    /* each n (values at the start) */
	nwint C;           /* number of change events */
	nwint Cp;          /* number of composition changes */
	double *pchg;      /* Cp x 3: actor (1..n), time, new presence (0/1) */
	double *pinit;     /* n values (presence at the start) or NULL */
	double *chg;       /* C x 4: covariate row (1..7; 15 = row 5 of a variable that only a windowed tertius uses), actor (1..n), time, new value */
	double *opp;       /* N x n or NULL */
	double *tie;       /* n x n or NULL */
	double *blob;
	int eff[64], P;    /* active effect slots (1-based; 100 + t = extra term t) and their number */
	int ign[20];       /* ignoreRep of the built-in effect slots */
	int nxl, nxt, nxa;
	XLayer xl[XMAXL];
	XTermSpec xt[XMAXT];
	const double *xattr[XMAXA];
} Spec;

typedef struct {
	double *V[2];                      /* the term's own copies of its networks (values) */
	nwint pc[2], pe[2];                /* next creation / expiry per network */
	double *C;                         /* cache (counts) or statistic (tie) */
	double *S;                         /* tertiusDiff: statistic matrix */
	double *cv;                        /* tertiusDiff: cache vector (NaN = no in-neighbour) */
	int impute;                        /* tertiusDiff: some cache is NaN */
} XTState;

typedef struct {
	nwint n;
	XTState *xs;
	int nxt, err;                      /* err: 1/2 = a situation in which goldfish's tertiusDiff stops with an error */
	double *tv;                        /* tie values, n x n */
	char *Ab;                          /* tv > 0 */
	double *inS, *outS, *inV, *outV;   /* column / row sums of sign(tv) and of tv */
	double *inA, *outA;                /* column / row counts of Ab */
	double *W[5];                      /* windowed tie values (slot 0 inertia, 1 recip, 2 indeg, 3 outdeg, 4 tertius) */
	double *tcache, *tstat;            /* windowed tertius: goldfish's cache (mean attribute of the in-neighbours, NaN = none) and statistic */
	double *wInS, *wInV, *wOutS, *wOutV;
	nwint ptr[5];
	double *cv[7];                     /* the covariates as they are now (changed by the change events) */
	double *cvbuf;
	nwint cpos;                        /* next change event */
	char *pres;                        /* presence of the actors now, NULL = nobody ever absent (no composition change) */
	nwint ppos;                        /* next composition change */
	int32_t *AA, *AtA, *AAt;
	double *nt, *tnum;
	int sym;
} State;

/* optional per-event output (code 30); NULL = a fit, which therefore does no extra work */
typedef struct {
	int write;                         /* 0 = only count the records */
	nwint nrec;                        /* records seen */
	nwint cap;                         /* rows allocated for the records */
	double *ll, *typ, *tim;            /* record rows */
	double *pr;                        /* probability blocks or NULL */
	nwint pstride;                     /* values per block */
} PredOut;
static PredOut *g_po = NULL;
static double g_lastS = 0.0;           /* rate model with intercept: sum of the hazards of the last exposure term */

static double sgn(double x) { return (double)((x > 0.0) - (x < 0.0)); }

static void state_free(State *st) {
	int q, t;
	if (st->xs) for (t = 0; t < st->nxt; t++) { free(st->xs[t].V[0]); free(st->xs[t].V[1]); free(st->xs[t].C); free(st->xs[t].S); free(st->xs[t].cv); }
	free(st->xs);
	free(st->tv); free(st->Ab); free(st->inS); free(st->outS); free(st->inV); free(st->outV);
	free(st->inA); free(st->outA);
	for (q = 0; q < 5; q++) free(st->W[q]);
	free(st->tcache); free(st->tstat);
	free(st->wInS); free(st->wInV); free(st->wOutS); free(st->wOutV);
	free(st->AA); free(st->AtA); free(st->AAt); free(st->nt); free(st->tnum); free(st->cvbuf); free(st->pres);
	memset(st, 0, sizeof(State));
}

static int xt_init(const Spec *sp, State *st);

static int has_eff(const Spec *sp, int e) {
	int k;
	for (k = 0; k < sp->P; k++) if (sp->eff[k] == e) return 1;
	return 0;
}

static int state_init(const Spec *sp, State *st) {
	nwint n = sp->n;
	size_t nn = (size_t)n * (size_t)n, un = (size_t)n;
	int q, ch = (sp->model == M_CHOICE), ra = (sp->model == M_RATE || sp->model == M_RATEINT), co = (sp->model == M_COORD);
	memset(st, 0, sizeof(State));
	st->n = n;
	st->sym = co;
	st->cvbuf = (double *)malloc(7 * un * sizeof(double));
	if (!st->cvbuf) return 909;
	for (q = 0; q < 7; q++) { st->cv[q] = st->cvbuf + (size_t)q * un; memcpy(st->cv[q], sp->cov[q], un * sizeof(double)); }
	st->tv = (double *)calloc(nn, sizeof(double));
	st->Ab = (char *)calloc(nn, 1);
	st->inS = (double *)calloc(un, sizeof(double)); st->outS = (double *)calloc(un, sizeof(double));
	st->inV = (double *)calloc(un, sizeof(double)); st->outV = (double *)calloc(un, sizeof(double));
	st->inA = (double *)calloc(un, sizeof(double)); st->outA = (double *)calloc(un, sizeof(double));
	if (!st->tv || !st->Ab || !st->inS || !st->outS || !st->inV || !st->outV || !st->inA || !st->outA) return 909;
	if (sp->Cp > 0 || sp->pinit) {
		nwint a;
		st->pres = (char *)malloc(un);
		if (!st->pres) return 909;
		for (a = 0; a < n; a++) st->pres[a] = (char)(sp->pinit ? (sp->pinit[a] != 0.0) : 1);
	}
	for (q = 0; q < 5; q++) {
		int use = (ch && q < 2) || (ra && q >= 2 && q < 4) || (co && (q == 0 || q == 2));
		if (q == 4) use = (ch && has_eff(sp, 16)) || (ra && has_eff(sp, 5)) || (co && has_eff(sp, 9));
		if (use && sp->win[q] < NOWIN) {
			st->W[q] = (double *)calloc(nn, sizeof(double));
			if (!st->W[q]) return 909;
		}
	}
	if (st->W[4]) {
		st->tcache = (double *)calloc(un, sizeof(double)); st->tstat = (double *)calloc(un, sizeof(double));
		if (!st->tcache || !st->tstat) return 909;
	}
	if ((ra || co) && sp->win[2] < NOWIN) { st->wInS = (double *)calloc(un, sizeof(double)); st->wInV = (double *)calloc(un, sizeof(double)); if (!st->wInS || !st->wInV) return 909; }
	if (ra && sp->win[3] < NOWIN) { st->wOutS = (double *)calloc(un, sizeof(double)); st->wOutV = (double *)calloc(un, sizeof(double)); if (!st->wOutS || !st->wOutV) return 909; }
	if (ch) {
		if (has_eff(sp, 9) || has_eff(sp, 10)) st->AA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 11) || has_eff(sp, 13)) st->AtA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 12)) st->AAt = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 14)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 16) && !st->W[4]) st->tnum = (double *)calloc(un, sizeof(double));
		if (((has_eff(sp, 9) || has_eff(sp, 10)) && !st->AA) || ((has_eff(sp, 11) || has_eff(sp, 13)) && !st->AtA) ||
			(has_eff(sp, 12) && !st->AAt) || (has_eff(sp, 14) && !st->nt) || (has_eff(sp, 16) && !st->W[4] && !st->tnum)) return 909;
	}
	else if (ra) {
		if (has_eff(sp, 4)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 5) && !st->W[4]) st->tnum = (double *)calloc(un, sizeof(double));
		if ((has_eff(sp, 4) && !st->nt) || (has_eff(sp, 5) && !st->W[4] && !st->tnum)) return 909;
	}
	else if (co) {
		/* coordination: the products are recomputed per event from Ab where needed; nodeTrans and tertius are kept */
		if (has_eff(sp, 7)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 9) && !st->W[4]) st->tnum = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 8)) st->AA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 10)) st->AtA = (int32_t *)calloc(nn, sizeof(int32_t));
		if ((has_eff(sp, 7) && !st->nt) || (has_eff(sp, 9) && !st->W[4] && !st->tnum) || (has_eff(sp, 8) && !st->AA) || (has_eff(sp, 10) && !st->AtA)) return 909;
	}
	return xt_init(sp, st);
}

/* a cell (a,b) of Ab changes by d = +1 / -1 (call BEFORE Ab is changed) */
static void flip_update(const Spec *sp, State *st, nwint a, nwint b, int d) {
	nwint n = st->n, i, j, k;
	const char *A = st->Ab;
	if (st->nt) {
		double s1 = 0.0, s2 = 0.0;
		for (k = 0; k < n; k++) {
			s1 += (double)(A[a * n + k] & A[k * n + b]);
			s2 += (double)(A[a * n + k] & A[b * n + k]);
		}
		st->nt[a] += d * (s1 + s2);
		for (i = 0; i < n; i++) if (A[i * n + a] && A[i * n + b]) st->nt[i] += d;
	}
	if (st->AA) {
		for (j = 0; j < n; j++) st->AA[a * n + j] += d * A[b * n + j];
		for (i = 0; i < n; i++) st->AA[i * n + b] += d * A[i * n + a];
	}
	if (st->AtA) {
		for (j = 0; j < n; j++) {
			if (j == b) continue;
			st->AtA[b * n + j] += d * A[a * n + j];
			st->AtA[j * n + b] += d * A[a * n + j];
		}
		st->AtA[b * n + b] += d;
	}
	if (st->AAt) {
		for (j = 0; j < n; j++) {
			if (j == a) continue;
			st->AAt[a * n + j] += d * A[j * n + b];
			st->AAt[j * n + a] += d * A[j * n + b];
		}
		st->AAt[a * n + a] += d;
	}
	st->inA[b] += d;
	st->outA[a] += d;
}

/* windowed tertius, a change of the tie (a,b) by dv: goldfish's update_REM_choice_tertius (cache of the mean attribute of the
   in-neighbours of every actor; NaN = none; the statistic of an actor without in-neighbour is imputed with the mean of the
   cache). Called BEFORE the tie value changes. */
static void tert_impute(State *st, int anynan_or_empty) {
	nwint j, n = st->n;
	double sum = 0.0, imp;
	int cnt = 0;
	if (!anynan_or_empty) return;
	for (j = 0; j < n; j++) if (!isnan(st->tcache[j])) { sum += st->tcache[j]; cnt++; }
	imp = cnt > 0 ? sum / (double)cnt : NAN;
	for (j = 0; j < n; j++) if (isnan(st->tcache[j])) st->tstat[j] = imp;
}

static void tert_tie(State *st, nwint a, nwint b, double dv) {
	nwint n = st->n, k, j;
	double old = st->W[4][a * n + b], nw = old + dv, sum = 0.0, val;
	int cnt = 0, anynan = 0, allzero = 1, so = (old > 0.0) - (old < 0.0), sn = (nw > 0.0) - (nw < 0.0);
	if (so == sn) return;
	for (j = 0; j < n; j++) { if (isnan(st->tcache[j])) anynan = 1; else if (st->tcache[j] != 0.0) allzero = 0; }
	for (k = 0; k < n; k++) if (k != a && st->W[4][k * n + b] > 0.0) { sum += st->cv[4][k]; cnt++; }
	if (sn - so == 1) { sum += st->cv[4][a]; cnt++; }
	val = cnt > 0 ? sum / (double)cnt : NAN;
	st->tcache[b] = val;
	if (!isnan(val)) st->tstat[b] = val;
	tert_impute(st, anynan || allzero || isnan(val));
}

/* windowed tertius, the attribute of actor a changes to nv (goldfish: the out-neighbours of a get a new mean) */
static void tert_attr(State *st, nwint a, double nv) {
	nwint n = st->n, k, j;
	int anynan = 0, allzero = 1;
	for (j = 0; j < n; j++) { if (isnan(st->tcache[j])) anynan = 1; else if (st->tcache[j] != 0.0) allzero = 0; }
	for (j = 0; j < n; j++) if (st->W[4][a * n + j] > 0.0) {
		double sum = nv;
		int cnt = 1;
		for (k = 0; k < n; k++) if (k != a && st->W[4][k * n + j] > 0.0) { sum += st->cv[4][k]; cnt++; }
		st->tcache[j] = sum / (double)cnt;
		st->tstat[j] = st->tcache[j];
	}
	tert_impute(st, anynan || allzero);
}

static void wupd(State *st, int q, nwint a, nwint b, double dv) {
	nwint n = st->n;
	double old = st->W[q][a * n + b], nw = old + dv;
	if (q == 4) tert_tie(st, a, b, dv);
	if (q == 2) { st->wInS[b] += sgn(nw) - sgn(old); st->wInV[b] += dv; }
	if (q == 3) { st->wOutS[a] += sgn(nw) - sgn(old); st->wOutV[a] += dv; }
	st->W[q][a * n + b] = nw;
}

static void apply_cell(const Spec *sp, State *st, nwint a, nwint b, double v) {
	nwint n = st->n;
	double old = st->tv[a * n + b], nw = old + v;
	int oa = (old > 0.0), na = (nw > 0.0), q;
	double so = sgn(old), sn = sgn(nw);
	st->inS[b] += sn - so; st->outS[a] += sn - so;
	st->inV[b] += v; st->outV[a] += v;
	if (oa != na) {
		flip_update(sp, st, a, b, na - oa);
		if (st->tnum) st->tnum[b] += (na - oa) * st->cv[4][a];
		st->Ab[a * n + b] = (char)na;
	}
	st->tv[a * n + b] = nw;
	for (q = 0; q < 5; q++) if (st->W[q]) wupd(st, q, a, b, v);
}

static void apply_event(const Spec *sp, State *st, nwint k) {
	nwint a = (nwint)sp->ev[5 * k] - 1, b = (nwint)sp->ev[5 * k + 1] - 1;
	double v = sp->ev[5 * k + 3];
	apply_cell(sp, st, a, b, v);
	if (st->sym) apply_cell(sp, st, b, a, v);
}

/* change events (time-varying covariates): apply the next change event */
static void apply_one_change(const Spec *sp, State *st) {
	const double *c = sp->chg + 4 * st->cpos;
	int q = ((int)c[0] % 10) - 1;
	nwint a = (nwint)c[1] - 1, b, n = st->n;
	double old = st->cv[q][a], nw = c[3];
	if (q == 4 && st->tnum) for (b = 0; b < n; b++) if (st->Ab[a * n + b]) st->tnum[b] += nw - old;
	if (q == 4 && st->W[4] && nw != old) tert_attr(st, a, nw);
	st->cv[q][a] = nw;
	st->cpos++;
}

/* all change events earlier than tcur (an event at tcur itself still sees the old value, as in goldfish) */
static void apply_changes(const Spec *sp, State *st, double tcur) {
	while (st->cpos < sp->C && sp->chg[4 * st->cpos + 2] < tcur) apply_one_change(sp, st);
}

/* composition changes: all changes at or before time T (goldfish: a change at the time of an event takes effect for that event) */
static void apply_pres(const Spec *sp, State *st, double T) {
	if (!st->pres) return;
	while (st->ppos < sp->Cp && sp->pchg[3 * st->ppos + 1] <= T) {
		st->pres[(nwint)sp->pchg[3 * st->ppos] - 1] = (char)(sp->pchg[3 * st->ppos + 2] != 0.0);
		st->ppos++;
	}
}

/* the earliest time at which an event of an earlier group (rows < gs) leaves a window, if it is before tcur (else 1e300) */
static double next_expiry(const Spec *sp, const State *st, nwint gs, double tcur) {
	int q;
	double b = 1e300;
	for (q = 0; q < 5; q++) {
		double e;
		if (!st->W[q] || st->ptr[q] >= gs) continue;
		e = sp->ev[5 * st->ptr[q] + 2] + sp->win[q];
		if (e < tcur && e < b) b = e;
	}
	return b;
}

/* drop the events of earlier groups (rows < gs) that have left their window at time b (event time + window <= b) */
static void expire_upto(const Spec *sp, State *st, nwint gs, double b) {
	int q;
	for (q = 0; q < 5; q++) {
		if (!st->W[q]) continue;
		while (st->ptr[q] < gs) {
			nwint p = st->ptr[q];
			if (sp->ev[5 * p + 2] + sp->win[q] > b) break;
			wupd(st, q, (nwint)sp->ev[5 * p] - 1, (nwint)sp->ev[5 * p + 1] - 1, -sp->ev[5 * p + 3]);
			if (st->sym) wupd(st, q, (nwint)sp->ev[5 * p + 1] - 1, (nwint)sp->ev[5 * p] - 1, -sp->ev[5 * p + 3]);
			st->ptr[q]++;
		}
	}
}

/* windows: drop the events of earlier groups (rows < gs) that are older than the window at time tcur */
static void expire(const Spec *sp, State *st, nwint gs, double tcur) {
	int q;
	nwint n = st->n;
	for (q = 0; q < 5; q++) {
		if (!st->W[q]) continue;
		while (st->ptr[q] < gs) {
			nwint p = st->ptr[q];
			if (sp->ev[5 * p + 2] + sp->win[q] >= tcur) break;   /* as goldfish: the expiry time is computed as time + window */
			wupd(st, q, (nwint)sp->ev[5 * p] - 1, (nwint)sp->ev[5 * p + 1] - 1, -sp->ev[5 * p + 3]);
			if (st->sym) wupd(st, q, (nwint)sp->ev[5 * p + 1] - 1, (nwint)sp->ev[5 * p] - 1, -sp->ev[5 * p + 3]);
			st->ptr[q]++;
		}
	}
	(void)n;
}


/* ------------------------------------------------------------------
   Extra terms: goldfish's update functions, translated literally (see the header)
   ------------------------------------------------------------------ */
static int xt_layer(const XTermSpec *ts, int m) { return m == 0 ? ts->la : ts->lb; }

static int xt_init(const Spec *sp, State *st) {
	nwint n = sp->n, i, j, k;
	size_t nn = (size_t)n * (size_t)n;
	int t, m;
	st->nxt = sp->nxt;
	if (sp->nxt == 0) return 0;
	st->xs = (XTState *)calloc((size_t)sp->nxt, sizeof(XTState));
	if (!st->xs) return 909;
	for (t = 0; t < sp->nxt; t++) {
		const XTermSpec *ts = &sp->xt[t];
		XTState *x = &st->xs[t];
		int two = (ts->lb >= 0), win = (ts->win < NOWIN), tdiff = (ts->type == XT_TDIFF);
		double *S1 = NULL, *S2 = NULL;
		for (m = 0; m < (two ? 2 : 1); m++) {
			int L = xt_layer(ts, m);
			x->V[m] = (double *)calloc(nn, sizeof(double));
			if (!x->V[m]) return 909;
			if (L > 0 && sp->xl[L - 1].hasmat) memcpy(x->V[m], sp->xl[L - 1].mat, nn * sizeof(double));
		}
		x->C = (double *)calloc(nn, sizeof(double));
		if (!x->C) return 909;
		if (tdiff || ts->type == XT_FOUR) {
			x->S = (double *)calloc(nn, sizeof(double));
			x->cv = (double *)calloc((size_t)n, sizeof(double));
			if (!x->S || !x->cv) return 909;
		}
		/* initial cache and statistic: from the start state of the network(s); windowed terms start at 0 (goldfish) */
		if (ts->type == XT_TIE) {
			if (!win) for (i = 0; i < (nwint)nn; i++) x->C[i] = ts->wgt != 0.0 ? x->V[0][i] : (double)(x->V[0][i] > 0.0);
			continue;
		}
		if (win) continue;
		S1 = (double *)malloc(nn * sizeof(double));
		if (!S1) return 909;
		for (i = 0; i < (nwint)nn; i++) S1[i] = sgn(x->V[0][i]);
		if (two) {
			S2 = (double *)malloc(nn * sizeof(double));
			if (!S2) { free(S1); return 909; }
			for (i = 0; i < (nwint)nn; i++) S2[i] = sgn(x->V[1][i]);
		}
		else S2 = S1;
		switch (ts->type) {
		case XT_TRANS: case XT_MTRANS: case XT_CYCLE: case XT_MCYCLE:
			/* P = S1 S2 ; cycle: transposed */
			for (i = 0; i < n; i++) for (k = 0; k < n; k++) {
				double a = S1[i * n + k];
				if (a == 0.0) continue;
				for (j = 0; j < n; j++) {
					double b = S2[k * n + j];
					if (b == 0.0) continue;
					if (ts->type == XT_TRANS || ts->type == XT_MTRANS) x->C[i * n + j] += a * b; else x->C[j * n + i] += a * b;
				}
			}
			if (ts->type == XT_CYCLE) for (i = 0; i < n; i++) x->C[i * n + i] = 0.0;
			break;
		case XT_CSEND: case XT_MCSEND:
			/* crossprod(S1, S2)[i,j] = sum_k S1[k,i] S2[k,j]; mixed: plus transpose */
			for (k = 0; k < n; k++) for (i = 0; i < n; i++) {
				double a = S1[k * n + i];
				if (a == 0.0) continue;
				for (j = 0; j < n; j++) {
					double b = S2[k * n + j];
					if (b == 0.0) continue;
					x->C[i * n + j] += a * b;
					if (ts->type == XT_MCSEND) x->C[j * n + i] += a * b;
				}
			}
			break;
		case XT_CRECV: case XT_MCRECV:
			/* tcrossprod(S1, S2)[i,j] = sum_k S1[i,k] S2[j,k]; mixed: plus transpose */
			for (i = 0; i < n; i++) for (k = 0; k < n; k++) {
				double a = S1[i * n + k];
				if (a == 0.0) continue;
				for (j = 0; j < n; j++) {
					double b = S2[j * n + k];
					if (b == 0.0) continue;
					x->C[i * n + j] += a * b;
					if (ts->type == XT_MCRECV) x->C[j * n + i] += a * b;
				}
			}
			break;
		case XT_FOUR: {
			/* stat = (S S') S - diag(rowSums(S*S)) S - S diag(colSums(S*S)) + S, S = sign(network) without its diagonal */
			double *T = (double *)calloc(nn, sizeof(double));
			int allz = 1;
			if (!T) { free(S1); return 909; }
			for (i = 0; i < n; i++) S1[i * n + i] = 0.0;
			for (i = 0; i < (nwint)nn; i++) if (S1[i] != 0.0) { allz = 0; break; }
			if (!allz) {
				for (i = 0; i < n; i++) for (k = 0; k < n; k++) {
					double a = S1[i * n + k];
					if (a == 0.0) continue;
					for (j = 0; j < n; j++) T[i * n + j] += a * S1[j * n + k];
				}
				for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
					double a = T[i * n + j];
					if (a == 0.0) continue;
					for (k = 0; k < n; k++) x->C[i * n + k] += a * S1[j * n + k];
				}
				/* rowSums / colSums of S*S: S is 0/+-1 so S*S = |S| */
				for (i = 0; i < n; i++) {
					double rs = 0.0;
					for (j = 0; j < n; j++) rs += S1[i * n + j] * S1[i * n + j];
					for (j = 0; j < n; j++) x->C[i * n + j] -= rs * S1[i * n + j];
				}
				for (j = 0; j < n; j++) {
					double cs = 0.0;
					for (i = 0; i < n; i++) cs += S1[i * n + j] * S1[i * n + j];
					for (i = 0; i < n; i++) x->C[i * n + j] -= S1[i * n + j] * cs;
				}
				for (i = 0; i < (nwint)nn; i++) x->C[i] += S1[i];
			}
			memcpy(x->S, x->C, nn * sizeof(double));
			free(T);
			break;
		}
		case XT_TDIFF: {
			const double *at = sp->xattr[ts->attr - 1];
			int allz = 1, any = 0;
			double sum = 0.0, imp;
			nwint cnt = 0;
			for (i = 0; i < (nwint)nn; i++) if (x->V[0][i] != 0.0) { allz = 0; break; }
			if (allz) break;       /* cache 0, statistic 0 */
			for (j = 0; j < n; j++) {
				double sa = 0.0; nwint c = 0;
				for (k = 0; k < n; k++) if (S1[k * n + j] == 1.0) { sa += at[k]; c++; }
				x->cv[j] = c > 0 ? sa / (double)c : NAN;
				if (c == 0) x->impute = 1;
			}
			for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
				double v = isnan(x->cv[j]) ? NAN : fabs(at[i] - x->cv[j]);
				x->S[i * n + j] = v;
				if (i != j && !isnan(v)) { sum += v; cnt++; any = 1; }
			}
			imp = any ? sum / (double)cnt : NAN;
			for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
				if (i == j) x->S[i * n + j] = 0.0;
				else if (isnan(x->S[i * n + j])) x->S[i * n + j] = imp;
			}
			break;
		}
		}
		free(S1);
		if (two) free(S2);
	}
	return 0;
}

/* tertiusDiff: the new in-neighbour mean of receiver r after the tie (s,r) changes by d = sign(new) - sign(old); called before V changes */
static void xt_tdiff(const Spec *sp, State *st, XTState *x, const XTermSpec *ts, nwint s, nwint r, int so, int sn) {
	nwint n = sp->n, i, j, k;
	const double *at = sp->xattr[ts->attr - 1], *V = x->V[0];
	int isEmpty = 1, isImpute, newval = sn - so, nchg = 0;
	double sum = 0.0, val, imp;
	nwint cnt = 0;
	{
		/* goldfish: isEmpty <- all(cache == 0), an NA when some cache is NaN and the others are all 0: "if (NA)" is an error */
		int anynan = 0, anynz = 0;
		for (j = 0; j < n; j++) { if (isnan(x->cv[j])) anynan = 1; else if (x->cv[j] != 0.0) anynz = 1; }
		if (anynan && !anynz) { st->err = 2; return; }
		isEmpty = !anynan && !anynz;
	}
	isImpute = 0;
	for (j = 0; j < n; j++) if (isnan(x->cv[j])) { isImpute = 1; break; }
	for (k = 0; k < n; k++) if (V[k * n + r] > 0.0 && !(newval != 1 && k == s)) { sum += at[k]; cnt++; }
	if (newval == 1) { sum += at[s]; cnt++; }
	val = cnt > 0 ? sum / (double)cnt : NAN;
	if (isnan(val) && ts->win < NOWIN) { st->err = 1; return; }   /* goldfish: with a window, a receiver without in-neighbour left makes its update fail ("subscript out of bounds") */
	if (!isImpute && isnan(val)) isImpute = 1;
	x->cv[r] = val;
	x->impute = isImpute;
	/* changes: column r for every actor but r */
	if (!isnan(val)) {
		for (i = 0; i < n; i++) if (i != r) x->S[i * n + r] = fabs(at[i] - val);
		nchg = 1;
	}
	if (isEmpty) {
		/* every cell that is not a change gets the mean of the changed values */
		if (nchg) {
			double ms = 0.0; nwint mc = 0;
			for (i = 0; i < n; i++) if (i != r) { ms += fabs(at[i] - val); mc++; }
			imp = mc > 0 ? ms / (double)mc : NAN;
		}
		else imp = NAN;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
			if (i == j) continue;
			if (nchg && j == r) continue;
			x->S[i * n + j] = imp;
		}
	}
	else if (isImpute) {
		/* statistic recomputed from the whole cache; its missing cells get its mean (diagonal excluded) */
		double ms = 0.0; nwint mc = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
			double v;
			if (i == j) continue;
			v = isnan(x->cv[j]) ? NAN : fabs(at[i] - x->cv[j]);
			if (!isnan(v)) { ms += v; mc++; }
		}
		imp = mc > 0 ? ms / (double)mc : NAN;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
			if (i == j) continue;
			if (isnan(x->cv[j]) && !(nchg && j == r)) x->S[i * n + j] = imp;
		}
	}
}

/* one change of network slot m of term t: cell (a,b) changes by dv (the term's own copy; goldfish: update_* with the network before the change) */
static void xt_cell(const Spec *sp, State *st, int t, int m, nwint a, nwint b, double dv) {
	nwint n = sp->n, k;
	const XTermSpec *ts = &sp->xt[t];
	XTState *x = &st->xs[t];
	double old = x->V[m][a * n + b], nw = old + dv, so = sgn(old), sn = sgn(nw), d = sn - so;
	const double *V0 = x->V[0], *V1 = x->V[ts->lb >= 0 ? 1 : 0];
	double *C = x->C;
	if (ts->type == XT_TIE) {
		x->V[m][a * n + b] = nw;
		if (ts->wgt != 0.0) C[a * n + b] = nw; else C[a * n + b] = (double)(nw > 0.0);
		return;
	}
	if (ts->type == XT_TDIFF) {
		if (a != b && d != 0.0) xt_tdiff(sp, st, x, ts, a, b, (int)so, (int)sn);
		x->V[m][a * n + b] = nw;
		return;
	}
	if (a == b || d == 0.0) { x->V[m][a * n + b] = nw; return; }
	if (ts->type == XT_FOUR) {
		/* goldfish update_DyNAM_choice_four: the paths are counted in the network after the change (as 0/1) */
		double *tmpv = (double *)malloc((size_t)n * sizeof(double)), *tmpw = (double *)malloc((size_t)n * sizeof(double));
		nwint i, j, l, pass;
		double inc = 2.0 * sn - 1.0;
		char *B = (char *)malloc((size_t)n * (size_t)n);
		if (!tmpv || !tmpw || !B) { free(tmpv); free(tmpw); free(B); return; }
		x->V[m][a * n + b] = nw;
		for (i = 0; i < n * n; i++) B[i] = (char)(V0[i] > 0.0);
		for (i = 0; i < n; i++) B[i * n + i] = 0;
		for (pass = 0; pass < 2; pass++) {
			/* pass 0 adds to the cache, pass 1 sets the statistic of the cells touched */
			for (i = 0; i < n; i++) {
				double v = 0.0, w = 0.0;
				for (k = 0; k < n; k++) {
					if (k != b && B[a * n + k] && B[i * n + k]) v += 1.0;
					if (k != a && B[k * n + b] && B[k * n + i]) w += 1.0;
				}
				tmpv[i] = (i == a) ? 0.0 : v;
				tmpw[i] = (i == b) ? 0.0 : w;
			}
			for (i = 0; i < n; i++) if (tmpv[i] > 0.0) { if (pass == 0) C[i * n + b] += inc * tmpv[i]; else x->S[i * n + b] = C[i * n + b] >= 0.0 ? C[i * n + b] : 0.0; }
			for (j = 0; j < n; j++) if (tmpw[j] > 0.0) { if (pass == 0) C[a * n + j] += inc * tmpw[j]; else x->S[a * n + j] = C[a * n + j] >= 0.0 ? C[a * n + j] : 0.0; }
			for (l = 0; l < n; l++) {
				if (l == a || !B[l * n + b]) continue;
				for (k = 0; k < n; k++) {
					if (k == b || !B[a * n + k]) continue;
					if (pass == 0) C[l * n + k] += inc; else x->S[l * n + k] = C[l * n + k] >= 0.0 ? C[l * n + k] : 0.0;
				}
			}
		}
		free(tmpv); free(tmpw); free(B);
		return;
	}
	switch (ts->type) {
	case XT_TRANS:
		for (k = 0; k < n; k++) {
			if (k != b && V0[b * n + k] > 0.0) C[a * n + k] += d;
			if (k != a && V0[k * n + a] > 0.0) C[k * n + b] += d;
		}
		break;
	case XT_CYCLE:
		for (k = 0; k < n; k++) {
			if (k == a || k == b) continue;
			if (V0[b * n + k] > 0.0) C[k * n + a] += d;
			if (V0[k * n + a] > 0.0) C[b * n + k] += d;
		}
		break;
	case XT_CSEND:
		for (k = 0; k < n; k++) if (k != a && k != b && V0[a * n + k] > 0.0) { C[k * n + b] += d; C[b * n + k] += d; }
		break;
	case XT_CRECV:
		for (k = 0; k < n; k++) if (k != a && k != b && V0[k * n + b] > 0.0) { C[a * n + k] += d; C[k * n + a] += d; }
		break;
	case XT_MTRANS:
		for (k = 0; k < n; k++) {
			if (k == a || k == b) continue;
			if (m == 0 ? V1[b * n + k] > 0.0 : V0[k * n + a] > 0.0) C[m == 0 ? a * n + k : k * n + b] += d;
		}
		break;
	case XT_MCYCLE:
		for (k = 0; k < n; k++) {
			if (k == a || k == b) continue;
			if (m == 0 ? V1[b * n + k] > 0.0 : V0[k * n + a] > 0.0) C[m == 0 ? k * n + a : b * n + k] += d;
		}
		break;
	case XT_MCSEND:
		for (k = 0; k < n; k++) {
			if (k == a || k == b) continue;
			if (m == 0 ? V1[a * n + k] > 0.0 : V0[a * n + k] > 0.0) { C[(m == 0 ? b * n + k : k * n + b)] += d; C[(m == 0 ? k * n + b : b * n + k)] += d; }
		}
		break;
	case XT_MCRECV:
		for (k = 0; k < n; k++) {
			if (k == a || k == b) continue;
			if (m == 0 ? V1[k * n + b] > 0.0 : V0[k * n + b] > 0.0) { C[(m == 0 ? a * n + k : k * n + a)] += d; C[(m == 0 ? k * n + a : a * n + k)] += d; }
		}
		break;
	}
	x->V[m][a * n + b] = nw;
}

static void xt_change(const Spec *sp, State *st, int t, int m, int sym, nwint a, nwint b, double dv) {
	xt_cell(sp, st, t, m, a, b, dv);
	if (sym) xt_cell(sp, st, t, m, b, a, dv);
}

/* bring every extra term up to the start of the group that begins at row gs of the dependent events (time tcur): the events of
   the networks with time < tcur (dependent events: rows < gs) and the expiries of windowed events with time + window < tcur,
   in the order of time; at equal times creations come before expiries, the first network before the second. */
static void xadvance(const Spec *sp, State *st, nwint gs, double tcur) {
	int t;
	for (t = 0; t < sp->nxt; t++) {
		const XTermSpec *ts = &sp->xt[t];
		XTState *x = &st->xs[t];
		int nsl = (ts->lb >= 0) ? 2 : 1, win = (ts->win < NOWIN);
		for (;;) {
			int bm = -1, bk = 0, m;
			double bt = 1e300;
			for (m = 0; m < nsl; m++) {
				int L = xt_layer(ts, m);
				nwint ne = (L == 0) ? sp->N : sp->xl[L - 1].ne;
				const double *ev = (L == 0) ? sp->ev : sp->xl[L - 1].ev;
				int stride = (L == 0) ? 5 : 4;
				if (x->pc[m] < ne) {
					double tc = ev[stride * x->pc[m] + 2];
					int ok = (L == 0) ? (x->pc[m] < gs) : (tc < tcur);
					if (ok && (bm < 0 || tc < bt)) { bm = m; bk = 0; bt = tc; }
				}
				if (win && x->pe[m] < x->pc[m]) {
					double te = ev[stride * x->pe[m] + 2] + ts->win;
					if (te < tcur && (bm < 0 || te < bt)) { bm = m; bk = 1; bt = te; }
				}
			}
			if (bm < 0) break;
			{
				int L = xt_layer(ts, bm);
				int stride = (L == 0) ? 5 : 4, sym = (L == 0) ? st->sym : sp->xl[L - 1].sym;
				const double *ev = (L == 0) ? sp->ev : sp->xl[L - 1].ev;
				nwint i = bk == 0 ? x->pc[bm] : x->pe[bm];
				nwint a = (nwint)ev[stride * i] - 1, b = (nwint)ev[stride * i + 1] - 1;
				double v = ev[stride * i + 3];
				xt_change(sp, st, t, bm, sym, a, b, bk == 0 ? v : -v);
				if (bk == 0) x->pc[bm]++; else x->pe[bm]++;
			}
		}
	}
}

static int eff_ign(const Spec *sp, int e) { return e >= 100 ? sp->xt[e - 100].ign : (e < 20 ? sp->ign[e] : 0); }

/* statistic of term t for the dyad (a,b) */
static double xt_stat(const Spec *sp, const State *st, int t, nwint a, nwint b) {
	const XTState *x = &st->xs[t];
	double v = (sp->xt[t].type == XT_TDIFF || sp->xt[t].type == XT_FOUR) ? x->S[a * sp->n + b] : x->C[a * sp->n + b];
	if (isnan(v)) v = 0.0;
	return v;
}

/* ------------------------------------------------------------------
   Statistics of one candidate vector. X is P x n (row k = k-th active effect)
   ------------------------------------------------------------------ */
static void choice_stats(const Spec *sp, const State *st, nwint s, double *X, double *tmp) {
	nwint n = sp->n, j, k, l;
	int q;
	const double *tv = st->tv;
	const char *A = st->Ab;
	for (q = 0; q < sp->P; q++) {
		double *x = X + (size_t)q * (size_t)n;
		switch (sp->eff[q]) {
		case 1:
			if (st->W[0]) {
				const double *w = st->W[0] + s * n;
				if (sp->wgt[0] != 0.0) for (j = 0; j < n; j++) x[j] = w[j]; else for (j = 0; j < n; j++) x[j] = (w[j] > 0.0);
			}
			else if (sp->wgt[0] != 0.0) for (j = 0; j < n; j++) x[j] = tv[s * n + j];
			else for (j = 0; j < n; j++) x[j] = (double)A[s * n + j];
			break;
		case 2:
			if (st->W[1]) {
				const double *w = st->W[1];
				if (sp->wgt[1] != 0.0) for (j = 0; j < n; j++) x[j] = w[j * n + s]; else for (j = 0; j < n; j++) x[j] = (w[j * n + s] > 0.0);
			}
			else if (sp->wgt[1] != 0.0) for (j = 0; j < n; j++) x[j] = tv[j * n + s];
			else for (j = 0; j < n; j++) x[j] = (double)A[j * n + s];
			break;
		case 3: for (j = 0; j < n; j++) x[j] = (sp->wgt[2] != 0.0) ? st->inV[j] : st->inS[j]; break;
		case 4: for (j = 0; j < n; j++) x[j] = (st->cv[0][j] == st->cv[0][s]); break;
		case 5: for (j = 0; j < n; j++) x[j] = fabs(st->cv[1][j] - st->cv[1][s]); break;
		case 6: for (j = 0; j < n; j++) x[j] = -fabs(st->cv[2][j] - st->cv[2][s]); break;
		case 7: for (j = 0; j < n; j++) x[j] = st->cv[3][j]; break;
		case 8: for (j = 0; j < n; j++) x[j] = (sp->wgt[3] != 0.0) ? st->outV[j] : st->outS[j]; break;
		case 9: for (j = 0; j < n; j++) x[j] = (double)st->AA[s * n + j]; break;
		case 10: for (j = 0; j < n; j++) x[j] = (double)st->AA[j * n + s]; break;
		case 11: for (j = 0; j < n; j++) x[j] = (double)st->AtA[s * n + j]; break;
		case 12: for (j = 0; j < n; j++) x[j] = (double)st->AAt[s * n + j]; break;
		case 13:
			for (j = 0; j < n; j++) tmp[j] = 0.0;
			for (k = 0; k < n; k++) if (A[s * n + k]) for (l = 0; l < n; l++) tmp[l] += (double)st->AtA[k * n + l];
			for (j = 0; j < n; j++) {
				double a = (double)A[s * n + j];
				x[j] = tmp[j] - a * (st->outA[s] + (double)st->AtA[j * n + j] - a);
			}
			break;
		case 14: for (j = 0; j < n; j++) x[j] = st->nt[j]; break;
		case 15: for (j = 0; j < n; j++) x[j] = st->cv[5][s] * st->cv[6][j]; break;
		case 16:
			if (st->W[4]) for (j = 0; j < n; j++) x[j] = isnan(st->tstat[j]) ? 0.0 : st->tstat[j];   /* goldfish turns a missing statistic into 0 */
			else for (j = 0; j < n; j++) x[j] = (st->inA[j] > 0.0) ? st->tnum[j] / st->inA[j] : 0.0;
			break;
		case 17: for (j = 0; j < n; j++) x[j] = sp->tie[s * n + j]; break;
		}
		if (sp->eff[q] >= 100) { int t = sp->eff[q] - 100; for (j = 0; j < n; j++) x[j] = xt_stat(sp, st, t, s, j); }
		if (eff_ign(sp, sp->eff[q])) for (j = 0; j < n; j++) if (A[s * n + j]) x[j] = 0.0;
	}
}

/* rate statistics (the same for every actor-event of a group): X is Pt x n, row 0 = ones when `icpt' */
static void rate_stats(const Spec *sp, const State *st, double *X, int icpt) {
	nwint n = sp->n, j;
	int q, off = icpt ? 1 : 0;
	if (icpt) for (j = 0; j < n; j++) X[j] = 1.0;
	for (q = 0; q < sp->P; q++) {
		double *x = X + (size_t)(q + off) * (size_t)n;
		switch (sp->eff[q]) {
		case 1:
			if (st->W[2]) for (j = 0; j < n; j++) x[j] = (sp->wgt[2] != 0.0) ? st->wInV[j] : st->wInS[j];
			else for (j = 0; j < n; j++) x[j] = (sp->wgt[2] != 0.0) ? st->inV[j] : st->inS[j];
			break;
		case 2:
			if (st->W[3]) for (j = 0; j < n; j++) x[j] = (sp->wgt[3] != 0.0) ? st->wOutV[j] : st->wOutS[j];
			else for (j = 0; j < n; j++) x[j] = (sp->wgt[3] != 0.0) ? st->outV[j] : st->outS[j];
			break;
		case 3: for (j = 0; j < n; j++) x[j] = st->cv[3][j]; break;
		case 4: for (j = 0; j < n; j++) x[j] = st->nt[j]; break;
		case 5:
			if (st->W[4]) for (j = 0; j < n; j++) x[j] = isnan(st->tstat[j]) ? 0.0 : st->tstat[j];   /* goldfish turns a missing statistic into 0 */
			else for (j = 0; j < n; j++) x[j] = (st->inA[j] > 0.0) ? st->tnum[j] / st->inA[j] : 0.0;
			break;
		}
	}
}

/* ------------------------------------------------------------------
   Evaluation: choice and rate (no intercept), softmax over a candidate set
   ------------------------------------------------------------------ */
static int eval_softmax(const Spec *sp, const double *theta, int todo, double *ll, double *g, double *H) {
	State st;
	nwint n = sp->n, N = sp->N, i, j, gs = 0;
	int P = sp->P, k, l, rc = 0, choice = (sp->model == M_CHOICE);
	double *X = NULL, *lp = NULL, *Pr = NULL, *m = NULL, *E = NULL, *tmp = NULL, mx, sumP, lrsum;
	double L = 0.0;
	int fresh = 1;

	if ((rc = state_init(sp, &st)) != 0) { state_free(&st); return rc; }
	X = (double *)malloc((size_t)P * (size_t)n * sizeof(double));
	lp = (double *)malloc((size_t)n * sizeof(double));
	Pr = (double *)malloc((size_t)n * sizeof(double));
	m = (double *)malloc((size_t)P * sizeof(double));
	E = (double *)malloc((size_t)P * (size_t)P * sizeof(double));
	tmp = (double *)malloc((size_t)n * sizeof(double));
	if (!X || !lp || !Pr || !m || !E || !tmp) { rc = 909; goto done; }
	for (k = 0; k < P; k++) g[k] = 0.0;
	if (todo >= 2) for (k = 0; k < P * P; k++) H[k] = 0.0;

	for (i = 0; i < N; i++) {
		nwint s = (nwint)sp->ev[5 * i] - 1, r = (nwint)sp->ev[5 * i + 1] - 1;
		if (fresh) {
			expire(sp, &st, gs, sp->ev[5 * i + 2]);
			apply_changes(sp, &st, sp->ev[5 * i + 2]);
			if (sp->nxt > 0) {
				xadvance(sp, &st, gs, sp->ev[5 * i + 2]);
				if (st.err) { SF_error("dynam_sim: tertiusdiff(): a candidate has no in-neighbour left in the network, where goldfish stops with an error\n"); rc = 198; goto done; }
			}
			if (!choice) rate_stats(sp, &st, X, 0);
		}
		if (st.pres) {
			apply_pres(sp, &st, sp->ev[5 * i + 2]);
			if (!st.pres[s] || !st.pres[r]) { SF_error("dynam_sim: an event has an absent sender or receiver\n"); rc = 198; goto done; }
		}
		if (choice) choice_stats(sp, &st, s, X, tmp);
		/* linear predictor */
		for (j = 0; j < n; j++) lp[j] = 0.0;
		for (k = 0; k < P; k++) {
			const double *x = X + (size_t)k * (size_t)n;
			double th = theta[k];
			for (j = 0; j < n; j++) lp[j] += th * x[j];
		}
		if (choice) {
			lp[s] = -1e300;
			if (sp->n1 > 0) for (j = 0; j < sp->n1; j++) lp[j] = -1e300;
			if (sp->opp) for (j = 0; j < n; j++) if (sp->opp[i * n + j] == 0.0) lp[j] = -1e300;
		}
		else if (sp->n1 > 0) for (j = sp->n1; j < n; j++) lp[j] = -1e300;
		if (st.pres) for (j = 0; j < n; j++) if (!st.pres[j]) lp[j] = -1e300;
		mx = lp[0];
		for (j = 1; j < n; j++) if (lp[j] > mx) mx = lp[j];
		sumP = 0.0;
		for (j = 0; j < n; j++) { Pr[j] = exp(lp[j] - mx); sumP += Pr[j]; }
		lrsum = mx + log(sumP);
		L += lp[choice ? r : s] - lrsum;
		for (j = 0; j < n; j++) Pr[j] /= sumP;
		if (g_po) {
			if (g_po->write && g_po->nrec < g_po->cap) {
				g_po->ll[g_po->nrec] = lp[choice ? r : s] - lrsum;
				g_po->typ[g_po->nrec] = (double)(i + 1);
				g_po->tim[g_po->nrec] = sp->ev[5 * i + 2];
				if (g_po->pr) for (j = 0; j < n; j++) g_po->pr[(size_t)i * (size_t)n + (size_t)j] = (lp[j] > -1e299) ? Pr[j] : SV_missval;
			}
			g_po->nrec++;
		}
		for (k = 0; k < P; k++) {
			const double *x = X + (size_t)k * (size_t)n;
			double a = 0.0;
			for (j = 0; j < n; j++) a += x[j] * Pr[j];
			m[k] = a;
			g[k] += x[choice ? r : s] - a;
		}
		if (todo >= 2) {
			for (k = 0; k < P * P; k++) E[k] = 0.0;
			for (j = 0; j < n; j++) {
				double pj = Pr[j];
				if (pj == 0.0) continue;
				for (k = 0; k < P; k++) {
					double xk = pj * X[(size_t)k * (size_t)n + j];
					for (l = k; l < P; l++) E[k * P + l] += xk * X[(size_t)l * (size_t)n + j];
				}
			}
			for (k = 0; k < P; k++) for (l = k; l < P; l++) {
				double h = -(E[k * P + l] - m[k] * m[l]);
				H[k * P + l] += h;
				if (l != k) H[l * P + k] += h;
			}
		}
		if (sp->ev[5 * i + 4] != 0.0) {
			nwint kk;
			for (kk = gs; kk <= i; kk++) apply_event(sp, &st, kk);
			gs = i + 1;
			fresh = 1;
		}
		else fresh = 0;
	}
	*ll = L;
done:
	free(X); free(lp); free(Pr); free(m); free(E); free(tmp);
	state_free(&st);
	return rc;
}

/* ------------------------------------------------------------------
   Evaluation: rate with intercept. theta[0] = intercept, then the P effects.
   ------------------------------------------------------------------ */
/* exposure term -w * sum_j exp(x_j theta) of a period of length w in the current state (X); 1 = overflow */
static int rateint_exposure(const Spec *sp, const double *theta, int todo, const double *X, double w, double *L, double *g, double *H,
		double *Sx, double *Sxx, const char *pres) {
	nwint n = sp->n, j, nrisk = (sp->n1 > 0 ? sp->n1 : sp->n);
	int Pt = sp->P + 1, k, l;
	double S = 0.0;
	{ const char *dp = getenv("DYNDBG"); if (dp && w > 0) { FILE *f = fopen(dp, "a"); if (f) { fprintf(f, "w %.6f :", w); for (j = 0; j < nrisk; j++) fprintf(f, " %.6f", X[(size_t)(Pt - 1) * (size_t)n + j]); fprintf(f, "\n"); fclose(f); } } }
	for (k = 0; k < Pt; k++) Sx[k] = 0.0;
	if (todo >= 2) for (k = 0; k < Pt * Pt; k++) Sxx[k] = 0.0;
	for (j = 0; j < nrisk; j++) {
		double lp = 0.0, ej;
		if (pres && !pres[j]) continue;
		for (k = 0; k < Pt; k++) lp += theta[k] * X[(size_t)k * (size_t)n + j];
		if (lp > 700.0) return 1;
		ej = exp(lp);
		S += ej;
		for (k = 0; k < Pt; k++) {
			double xk = ej * X[(size_t)k * (size_t)n + j];
			Sx[k] += xk;
			if (todo >= 2) for (l = k; l < Pt; l++) Sxx[k * Pt + l] += xk * X[(size_t)l * (size_t)n + j];
		}
	}
	g_lastS = S;
	*L -= w * S;
	for (k = 0; k < Pt; k++) g[k] -= w * Sx[k];
	if (todo >= 2) for (k = 0; k < Pt; k++) for (l = k; l < Pt; l++) {
		double h = -w * Sxx[k * Pt + l];
		H[k * Pt + l] += h;
		if (l != k) H[l * Pt + k] += h;
	}
	return 0;
}

/* a right-censored interval of the rate model with intercept (goldfish lists it as an event of its own) */
static void cens_record(double t, double w) {
	if (!g_po) return;
	if (g_po->write && g_po->nrec < g_po->cap) {
		g_po->ll[g_po->nrec] = -w * g_lastS;
		g_po->typ[g_po->nrec] = 0.0;
		g_po->tim[g_po->nrec] = t;
	}
	g_po->nrec++;
}

static int eval_rateint(const Spec *sp, const double *theta, int todo, double *ll, double *g, double *H, int *overflow) {
	State st;
	nwint n = sp->n, N = sp->N, i, gs = 0;
	int Pt = sp->P + 1, k, rc = 0;
	double *X = NULL, *Sx = NULL, *Sxx = NULL, L = 0.0, tprev = sp->ev[2], tlastchg = -1e300;
	int fresh = 1;

	*overflow = 0;
	/* goldfish observes from the earliest event of any list: a change before the first event starts the observation */
	{
		/* (changes of a variable that only a windowed tertius uses, code 15, do not set the start or the end, as in goldfish) */
		nwint ci;
		for (ci = 0; ci < sp->C; ci++) if (sp->chg[4 * ci] < 10.0) { if (sp->chg[4 * ci + 2] < tprev) tprev = sp->chg[4 * ci + 2]; break; }
		for (ci = sp->C - 1; ci >= 0; ci--) if (sp->chg[4 * ci] < 10.0) { tlastchg = sp->chg[4 * ci + 2]; break; }
	}
	if ((rc = state_init(sp, &st)) != 0) { state_free(&st); return rc; }
	X = (double *)malloc((size_t)Pt * (size_t)n * sizeof(double));
	Sx = (double *)malloc((size_t)Pt * sizeof(double));
	Sxx = (double *)malloc((size_t)Pt * (size_t)Pt * sizeof(double));
	if (!X || !Sx || !Sxx) { rc = 909; goto done; }
	for (k = 0; k < Pt; k++) g[k] = 0.0;
	if (todo >= 2) for (k = 0; k < Pt * Pt; k++) H[k] = 0.0;
	/* the hazard sums (S, Sx, Sxx) depend on the state only: computed once per period of constant state; a group of tied events has
	   the period before it (w > 0) at its first event only. Change events inside that period split it (the covariates of the
	   risk set change there). */
	for (i = 0; i < N; i++) {
		nwint s = (nwint)sp->ev[5 * i] - 1;
		double tcur = sp->ev[5 * i + 2], lps = 0.0, w;
		if (fresh) {
			/* no expire() here: events that leave a window inside the period are handled at their expiry time below */
			rate_stats(sp, &st, X, 1);
			{
				/* periods ended by a change event earlier than this event (before the first event: the observation
				   starts with the first change, as goldfish) */
				double t0 = tprev;
				for (;;) {
					double tc = (st.cpos < sp->C && sp->chg[4 * st.cpos + 2] < tcur) ? sp->chg[4 * st.cpos + 2] : 1e300;
					double te = next_expiry(sp, &st, gs, tcur), tb = (tc < te ? tc : te);
					if (tb >= 1e300) break;
					if (tb > t0) {
						apply_pres(sp, &st, tb);
						if (rateint_exposure(sp, theta, todo, X, tb - t0, &L, g, H, Sx, Sxx, st.pres)) { *overflow = 1; goto done; }
						cens_record(tb, tb - t0);
						t0 = tb;
					}
					expire_upto(sp, &st, gs, tb);
					while (st.cpos < sp->C && sp->chg[4 * st.cpos + 2] == tb) apply_one_change(sp, &st);
					rate_stats(sp, &st, X, 1);
				}
				tprev = t0;
			}
		}
		w = tcur - tprev;
		if (st.pres) {
			apply_pres(sp, &st, tcur);
			if (!st.pres[s]) { SF_error("dynam_sim: an event has an absent sender\n"); rc = 198; goto done; }
		}
		for (k = 0; k < Pt; k++) lps += theta[k] * X[(size_t)k * (size_t)n + s];
		if (lps > 700.0) { *overflow = 1; goto done; }
		L += lps;
		for (k = 0; k < Pt; k++) g[k] += X[(size_t)k * (size_t)n + s];
		g_lastS = 0.0;
		if (w != 0.0) {
			if (rateint_exposure(sp, theta, todo, X, w, &L, g, H, Sx, Sxx, st.pres)) { *overflow = 1; goto done; }
		}
		if (g_po) {
			if (g_po->write && g_po->nrec < g_po->cap) {
				nwint jj, nrisk = (sp->n1 > 0 ? sp->n1 : sp->n);
				g_po->ll[g_po->nrec] = lps - w * g_lastS;
				g_po->typ[g_po->nrec] = (double)(i + 1);
				g_po->tim[g_po->nrec] = tcur;
				if (g_po->pr) for (jj = 0; jj < n; jj++) {
					double lpj = 0.0;
					if (jj >= nrisk || (st.pres && !st.pres[jj])) { g_po->pr[(size_t)i * (size_t)n + (size_t)jj] = SV_missval; continue; }
					for (k = 0; k < Pt; k++) lpj += theta[k] * X[(size_t)k * (size_t)n + jj];
					g_po->pr[(size_t)i * (size_t)n + (size_t)jj] = lpj - w * g_lastS;
				}
			}
			g_po->nrec++;
		}
		tprev = tcur;
		if (sp->ev[5 * i + 4] != 0.0) {
			nwint kk;
			for (kk = gs; kk <= i; kk++) apply_event(sp, &st, kk);
			gs = i + 1;
			fresh = 1;
		}
		else fresh = 0;
	}
	/* ... and the last change after the last event ends it: the periods up to each remaining change count (and up to the
	   expiry of events from their windows before that change) */
	if (st.cpos < sp->C && tlastchg > tprev) {
		double t0 = tprev, tend = tlastchg;
		rate_stats(sp, &st, X, 1);
		while (st.cpos < sp->C && sp->chg[4 * st.cpos + 2] <= tend) {
			double tc = sp->chg[4 * st.cpos + 2];
			double te = next_expiry(sp, &st, gs, tend), tb = (tc < te ? tc : te);
			if (tb > t0) {
				apply_pres(sp, &st, tb);
				if (rateint_exposure(sp, theta, todo, X, tb - t0, &L, g, H, Sx, Sxx, st.pres)) { *overflow = 1; goto done; }
				cens_record(tb, tb - t0);
				t0 = tb;
			}
			expire_upto(sp, &st, gs, tb);
			while (st.cpos < sp->C && sp->chg[4 * st.cpos + 2] == tb) apply_one_change(sp, &st);
			rate_stats(sp, &st, X, 1);
		}
	}
	*ll = L;
done:
	free(X); free(Sx); free(Sxx);
	state_free(&st);
	return rc;
}


/* ------------------------------------------------------------------
   Evaluation: choice_coordination (undirected ties). Effect slots: 1 inertia, 2 indeg, 3 same, 4 diff, 5 sim, 6 alter,
   7 nodeTrans, 8 trans, 9 tertius, 10 four, 11 egoAlterInt, 12 tie.
   For every event the n x n matrix of linear predictors, p(a,b) = row-normalised softmax without the diagonal,
   Pnorm(a,b) = p(a,b) p(b,a) / (sum / 2); log likelihood ln Pnorm(s,r); gradient from the row scores (see the reference
   evaluator). Gradient only (no Hessian).
   ------------------------------------------------------------------ */
static double coord_S0(const Spec *sp, const State *st, int e, nwint a, nwint b, const int32_t *AAtA);
static double coord_S(const Spec *sp, const State *st, int e, nwint a, nwint b, const int32_t *AAtA) {
	double v;
	if (e >= 100) v = xt_stat(sp, st, e - 100, a, b);
	else v = coord_S0(sp, st, e, a, b, AAtA);
	if (eff_ign(sp, e) && st->Ab[a * sp->n + b]) v = 0.0;
	return v;
}
static double coord_S0(const Spec *sp, const State *st, int e, nwint a, nwint b, const int32_t *AAtA) {
	nwint n = sp->n;
	switch (e) {
	case 1:
		if (st->W[0]) return (sp->wgt[0] != 0.0) ? st->W[0][a * n + b] : (double)(st->W[0][a * n + b] > 0.0);
		return (sp->wgt[0] != 0.0) ? st->tv[a * n + b] : (double)st->Ab[a * n + b];
	case 2:
		if (st->W[2]) return (sp->wgt[2] != 0.0) ? st->wInV[b] : st->wInS[b];
		return (sp->wgt[2] != 0.0) ? st->inV[b] : st->inS[b];
	case 3: return (double)(st->cv[0][a] == st->cv[0][b]);
	case 4: return fabs(st->cv[1][a] - st->cv[1][b]);
	case 5: return -fabs(st->cv[2][a] - st->cv[2][b]);
	case 6: return st->cv[3][b];
	case 7: return st->nt[b];
	case 8: return (double)st->AA[a * n + b];
	case 9:
		if (st->W[4]) return isnan(st->tstat[b]) ? 0.0 : st->tstat[b];
		return (st->inA[b] > 0.0) ? st->tnum[b] / st->inA[b] : 0.0;
	case 10: {
		double A = (double)st->Ab[a * n + b];
		return (double)AAtA[a * n + b] - A * (st->outA[a] + (double)st->AtA[b * n + b] - A);
	}
	case 11: return st->cv[5][a] * st->cv[6][b];
	case 12: return sp->tie[a * n + b];
	}
	return 0.0;
}

static int eval_coord(const Spec *sp, const double *theta, int todo, double *ll, double *g) {
	State st;
	nwint n = sp->n, N = sp->N, i, a, b, k, l, gs = 0;
	int P = sp->P, q, rc = 0, dirty = 1, four = has_eff(sp, 10);
	double *stm = NULL, *p = NULL, *Pn = NULL, *rowexp = NULL, *rowPn = NULL, L = 0.0;
	int32_t *AAtA = NULL;
	(void)todo;

	if ((rc = state_init(sp, &st)) != 0) { state_free(&st); return rc; }
	stm = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	p = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	Pn = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
	rowexp = (double *)malloc((size_t)n * sizeof(double));
	rowPn = (double *)malloc((size_t)n * sizeof(double));
	if (four) AAtA = (int32_t *)calloc((size_t)n * (size_t)n, sizeof(int32_t));
	if (!stm || !p || !Pn || !rowexp || !rowPn || (four && !AAtA)) { rc = 909; goto done; }
	for (q = 0; q < P; q++) g[q] = 0.0;

	for (i = 0; i < N; i++) {
		nwint s = (nwint)sp->ev[5 * i] - 1, r = (nwint)sp->ev[5 * i + 1] - 1;
		double tot, Z;
		if (dirty) {
			expire(sp, &st, gs, sp->ev[5 * i + 2]);
			apply_changes(sp, &st, sp->ev[5 * i + 2]);
			if (sp->nxt > 0) {
				xadvance(sp, &st, gs, sp->ev[5 * i + 2]);
				if (st.err) { SF_error("dynam_sim: tertiusdiff(): a candidate has no in-neighbour left in the network, where goldfish stops with an error\n"); rc = 198; goto done; }
			}
			if (four) {
				/* AAtA = A * (A'A), recomputed when the tie matrix changed */
				memset(AAtA, 0, (size_t)n * (size_t)n * sizeof(int32_t));
				for (a = 0; a < n; a++) for (k = 0; k < n; k++) if (st.Ab[a * n + k]) {
					const int32_t *row = st.AtA + k * n;
					int32_t *out = AAtA + a * n;
					for (b = 0; b < n; b++) out[b] += row[b];
				}
			}
			dirty = 0;
		}
		if (st.pres) {
			apply_pres(sp, &st, sp->ev[5 * i + 2]);
			if (!st.pres[s] || !st.pres[r]) { SF_error("dynam_sim: an event has an absent actor\n"); rc = 198; goto done; }
		}
		for (a = 0; a < n; a++) for (b = 0; b < n; b++) {
			double v = 0.0;
			for (q = 0; q < P; q++) v += theta[q] * coord_S(sp, &st, sp->eff[q], a, b, AAtA);
			stm[a * n + b] = v;
		}
		for (a = 0; a < n; a++) {
			double rmax = stm[a * n], rs = 0.0;
			for (b = 1; b < n; b++) if (stm[a * n + b] > rmax) rmax = stm[a * n + b];
			for (b = 0; b < n; b++) {
				double e = (a == b) ? 0.0 : exp(stm[a * n + b] - rmax);
				if (st.pres && (!st.pres[a] || !st.pres[b])) e = 0.0;
				p[a * n + b] = e; rs += e;
			}
			if (rs > 0.0) for (b = 0; b < n; b++) p[a * n + b] /= rs;
		}
		tot = 0.0;
		for (a = 0; a < n; a++) for (b = 0; b < n; b++) { double v = p[a * n + b] * p[b * n + a]; Pn[a * n + b] = v; tot += v; }
		Z = tot / 2.0;
		for (a = 0; a < n; a++) { double rs = 0.0; for (b = 0; b < n; b++) { Pn[a * n + b] /= Z; rs += Pn[a * n + b]; } rowPn[a] = rs; }
		L += log(Pn[s * n + r]);
		if (g_po) {
			if (g_po->write && g_po->nrec < g_po->cap) {
				g_po->ll[g_po->nrec] = log(Pn[s * n + r]);
				g_po->typ[g_po->nrec] = (double)(i + 1);
				g_po->tim[g_po->nrec] = sp->ev[5 * i + 2];
				if (g_po->pr) {
					size_t pos = (size_t)i * (size_t)g_po->pstride;
					for (a = 0; a < n; a++) for (b = a + 1; b < n; b++)
						g_po->pr[pos++] = (!st.pres || (st.pres[a] && st.pres[b])) ? Pn[a * n + b] : SV_missval;
				}
			}
			g_po->nrec++;
		}
		for (q = 0; q < P; q++) {
			int e = sp->eff[q];
			double T1 = 0.0, T2 = 0.0;
			for (a = 0; a < n; a++) {
				double re = 0.0;
				for (b = 0; b < n; b++) {
					double sv = coord_S(sp, &st, e, a, b, AAtA);
					re += p[a * n + b] * sv;
					T1 += Pn[a * n + b] * sv;
				}
				rowexp[a] = re;
				T2 += re * rowPn[a];
			}
			{
				double Ssr = coord_S(sp, &st, e, s, r, AAtA), Srs = coord_S(sp, &st, e, r, s, AAtA);
				g[q] += (Ssr - rowexp[s]) + (Srs - rowexp[r]) - (T1 - T2);
			}
		}
		if (sp->ev[5 * i + 4] != 0.0) {
			nwint kk;
			for (kk = gs; kk <= i; kk++) apply_event(sp, &st, kk);
			gs = i + 1;
			dirty = 1;
		}
		(void)l;
	}
	*ll = L;
done:
	free(stm); free(p); free(Pn); free(rowexp); free(rowPn); free(AAtA);
	state_free(&st);
	return rc;
}

/* ------------------------------------------------------------------
   Blob reading
   ------------------------------------------------------------------ */
static int spec_read(Spec *sp, nwint L) {
	nwint i, pos;
	int k;
	memset(sp, 0, sizeof(Spec));
	sp->blob = (double *)malloc((size_t)L * sizeof(double));
	if (!sp->blob) return 909;
	for (i = 0; i < L; i++) {
		ST_double v;
		SF_vdata(1, (ST_int)(i + 1), &v);
		sp->blob[i] = v;
	}
	if (L < HDR || (int)sp->blob[0] != DYNAM_PROTOCOL) return 198;
	sp->model = (int)sp->blob[1];
	sp->n = (nwint)sp->blob[2];
	sp->N = (nwint)sp->blob[3];
	sp->n1 = (nwint)sp->blob[4];
	sp->M = (int)sp->blob[5];
	sp->hasopp = (int)sp->blob[6];
	sp->hastie = (int)sp->blob[7];
	for (k = 0; k < 4; k++) { sp->win[k] = sp->blob[8 + k]; sp->wgt[k] = sp->blob[12 + k]; }
	sp->win[4] = sp->blob[39];
	for (k = 0; k < 20; k++) sp->act[k] = sp->blob[16 + k];
	sp->P = 0;
	for (k = 0; k < sp->M && k < 20; k++) if (sp->act[k] != 0.0) { sp->eff[sp->P++] = k + 1; sp->ign[k + 1 < 20 ? k + 1 : 19] = (sp->act[k] == 2.0); }
	pos = HDR;
	sp->C = (nwint)sp->blob[36];
	if (sp->C < 0) return 198;
	sp->ev = sp->blob + pos; pos += 5 * sp->N;
	sp->chg = sp->blob + pos; pos += 4 * sp->C;
	sp->Cp = (nwint)sp->blob[37];
	if (sp->Cp < 0) return 198;
	sp->pchg = sp->blob + pos; pos += 3 * sp->Cp;
	if (sp->blob[38] != 0.0) { sp->pinit = sp->blob + pos; pos += sp->n; }
	for (k = 0; k < 7; k++) { sp->cov[k] = sp->blob + pos; pos += sp->n; }
	if (sp->hasopp) { sp->opp = sp->blob + pos; pos += sp->N * sp->n; }
	if (sp->hastie) { sp->tie = sp->blob + pos; pos += sp->n * sp->n; }
	{
		nwint xlen = (nwint)sp->blob[40];
		if (xlen < 0) return 198;
		if (xlen > 0) {
			const double *xb = sp->blob + pos;
			nwint p = 3, nl, nt, na, l, t;
			if (xlen < 3) return 198;
			nl = (nwint)xb[0]; nt = (nwint)xb[1]; na = (nwint)xb[2];
			if (nl < 0 || nl > XMAXL || nt < 0 || nt > XMAXT || na < 0 || na > XMAXA) { SF_error("dynam_sim: too many exogenous networks or terms\n"); return 198; }
			for (l = 0; l < nl; l++) {
				XLayer *xl = &sp->xl[l];
				if (p + 3 > xlen) return 198;
				xl->hasmat = (int)xb[p]; xl->sym = (int)xb[p + 1]; xl->ne = (nwint)xb[p + 2]; p += 3;
				if (xl->ne < 0) return 198;
				xl->ev = xb + p; p += 4 * xl->ne;
				if (xl->hasmat) { xl->mat = xb + p; p += sp->n * sp->n; }
				if (p > xlen) return 198;
			}
			for (t = 0; t < nt; t++) {
				XTermSpec *ts = &sp->xt[t];
				if (p + 8 > xlen) return 198;
				ts->type = (int)xb[p]; ts->la = (int)xb[p + 1]; ts->lb = (int)xb[p + 2]; ts->win = xb[p + 3]; ts->wgt = xb[p + 4];
				ts->ign = (int)xb[p + 5]; ts->attr = (int)xb[p + 6]; p += 8;
				if (ts->type < XT_TIE || ts->type > XT_FOUR || ts->la < 0 || ts->la > nl || ts->lb > nl || ts->lb < -1) return 198;
				if (ts->type == XT_TDIFF && (ts->attr < 1 || ts->attr > na)) return 198;
				if (sp->P >= 60) return 198;
				sp->eff[sp->P++] = 100 + (int)t;
			}
			for (l = 0; l < na; l++) { if (p + sp->n > xlen) return 198; sp->xattr[l] = xb + p; p += sp->n; }
			if (p != xlen) return 198;
			sp->nxl = (int)nl; sp->nxt = (int)nt; sp->nxa = (int)na;
			pos += xlen;
		}
	}
	if (pos != L) return 198;
	if (sp->n < 2 || sp->N < 1) return 198;
	for (i = 0; i < sp->C; i++) {
		double r = sp->chg[4 * i], a = sp->chg[4 * i + 1];
		if (!((r >= 1.0 && r <= 7.0) || r == 15.0) || !( a >= 1.0 && a <= (double)sp->n)) {
			SF_error("dynam_sim: a change event has an actor or covariate outside the range\n");
			return 198;
		}
		if (i > 0 && sp->chg[4 * i + 2] < sp->chg[4 * i - 2]) { SF_error("dynam_sim: change events are not sorted by time\n"); return 198; }
	}
	for (i = 0; i < sp->Cp; i++) {
		double a = sp->pchg[3 * i], v = sp->pchg[3 * i + 2];
		if (!(a >= 1.0 && a <= (double)sp->n) || !(v == 0.0 || v == 1.0)) {
			SF_error("dynam_sim: a composition change has an actor outside 1..n or a presence other than 0/1\n");
			return 198;
		}
		if (i > 0 && sp->pchg[3 * i + 1] < sp->pchg[3 * i - 2]) { SF_error("dynam_sim: composition changes are not sorted by time\n"); return 198; }
	}
	if (sp->nxt > 0 && (sp->model == M_RATE || sp->model == M_RATEINT)) { SF_error("dynam_sim: extra terms are not available in the rate model\n"); return 198; }
	{
		int l;
		for (l = 0; l < sp->nxl; l++) {
			nwint e;
			for (e = 0; e < sp->xl[l].ne; e++) {
				double a = sp->xl[l].ev[4 * e], b = sp->xl[l].ev[4 * e + 1];
				if (!(a >= 1.0 && a <= (double)sp->n && b >= 1.0 && b <= (double)sp->n)) { SF_error("dynam_sim: an event of an exogenous network has an actor outside 1..n\n"); return 198; }
				if (e > 0 && sp->xl[l].ev[4 * e + 2] < sp->xl[l].ev[4 * e - 2]) { SF_error("dynam_sim: events of an exogenous network are not sorted by time\n"); return 198; }
			}
		}
	}
	/* every actor index of an event must be an actor 1..n (the kernels index arrays with it) */
	for (i = 0; i < sp->N; i++) {
		double a = sp->ev[5 * i], b = sp->ev[5 * i + 1];
		if (!(a >= 1.0 && a <= (double)sp->n && b >= 1.0 && b <= (double)sp->n)) {
			SF_error("dynam_sim: an event has a sender or receiver outside 1..n\n");
			return 198;
		}
	}
	return 0;
}

/* ------------------------------------------------------------------
   Series tools (code 40): the outlier rules of goldfish's examineOutliers (Hampel, IQR, Top) and the change point methods
   of goldfish's examineChangepoints, which calls cpt.mean() / cpt.var() of the R package 'changepoint' (Killick, Haynes,
   Eckley, Fearnhead, Lee; GPL, usable under GPL-3). The change point code below is a port of that package's C and R
   sources (version 2.3): PELT (PELT_one_func_minseglen.c), binary segmentation (BinSeg_one_func_minseglen.c),
   at most one change (single.norm.R), the cost functions (cost_general_functions.c) and the penalty (penalty_decision.R),
   all with the package's default penalty MBIC and test statistic Normal. Only what examineChangepoints can reach is ported.
   ------------------------------------------------------------------ */
#define SER_PI 3.141592653589793238462643383279502884

static int dcmp(const void *a, const void *b) {
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

/* R's median() of a vector (the buffer is sorted) */
static double med_sorted(double *v, nwint m) {
	qsort(v, (size_t)m, sizeof(double), dcmp);
	if (m % 2 == 1) return v[m / 2];
	return (v[m / 2 - 1] + v[m / 2]) / 2.0;
}

/* R's quantile(type = 7) of a SORTED vector */
static double quant7(const double *v, nwint m, double p) {
	double index = 1.0 + (double)(m > 1 ? m - 1 : 0) * p, lo = floor(index), hi = ceil(index), qs;
	qs = v[(nwint)lo - 1];
	if (index > lo && v[(nwint)hi - 1] != qs) {
		double h = index - lo;
		qs = (1.0 - h) * qs + h * v[(nwint)hi - 1];
	}
	return qs;
}

typedef double (*CptCost)(double, double, double, int);
static double cost_mean_mbic(double x, double x2, double x3, int n) { (void)x3; return x2 - (x * x) / n + log((double)n); }
static double cost_var_mbic(double x, double x2, double x3, int n) {
	(void)x; (void)x2;
	if (x3 <= 0) x3 = 0.00000000001;
	return n * (log(2 * SER_PI) + log(x3 / n) + 1) + log((double)n);
}

/* PELT (changepoint's PELTC): returns the change points (including n) sorted in cp, their number in *ncp */
static int cpt_pelt(const double *s1, const double *s2, const double *s3, int n, double pen, int minseglen, CptCost cf, int *cp, int *ncp) {
	double *lastchangelike = (double *)calloc((size_t)n + 1, sizeof(double)), *tmplike = (double *)calloc((size_t)n + 1, sizeof(double));
	int *lastchangecpts = (int *)calloc((size_t)n + 1, sizeof(int)), *checklist = (int *)calloc((size_t)n + 1, sizeof(int));
	int nchecklist, tstar, i, j, whichout, nchecktmp, last, k = 0;
	double minout;
	if (!lastchangelike || !tmplike || !lastchangecpts || !checklist) { free(lastchangelike); free(tmplike); free(lastchangecpts); free(checklist); return 909; }
	lastchangelike[0] = -pen;
	for (j = minseglen; j < 2 * minseglen; j++) lastchangelike[j] = cf(s1[j], s2[j], s3[j], j);
	nchecklist = 2;
	checklist[0] = 0;
	checklist[1] = minseglen;
	for (tstar = 2 * minseglen; tstar < n + 1; tstar++) {
		if (lastchangelike[tstar] == 0) {
			for (i = 0; i < nchecklist; i++)
				tmplike[i] = lastchangelike[checklist[i]] + cf(s1[tstar] - s1[checklist[i]], s2[tstar] - s2[checklist[i]], s3[tstar] - s3[checklist[i]], tstar - checklist[i]) + pen;
			minout = tmplike[0]; whichout = 0;
			for (i = 1; i < nchecklist; i++) if (tmplike[i] < minout) { minout = tmplike[i]; whichout = i; }
			lastchangelike[tstar] = minout;
			lastchangecpts[tstar] = checklist[whichout];
			nchecktmp = 0;
			for (i = 0; i < nchecklist; i++) if (tmplike[i] <= (lastchangelike[tstar] + pen)) checklist[nchecktmp++] = checklist[i];
			nchecklist = nchecktmp;
		}
		checklist[nchecklist] = tstar - (minseglen - 1);
		nchecklist++;
	}
	last = n;
	while (last != 0) { cp[k++] = last; last = lastchangecpts[last]; }
	/* ascending */
	for (i = 0; i < k / 2; i++) { int t = cp[i]; cp[i] = cp[k - 1 - i]; cp[k - 1 - i] = t; }
	*ncp = k;
	free(lastchangelike); free(tmplike); free(lastchangecpts); free(checklist);
	return 0;
}

/* binary segmentation (changepoint's binseg, Q = maximal number of change points); change points (including n) in cp */
static int cpt_binseg(const double *s1, const double *s2, const double *s3, int n, double pen, int minseglen, int Q, CptCost cf, int *cp, int *ncp) {
	double *lambda = (double *)malloc(((size_t)n + 1) * sizeof(double)), *likeout = (double *)calloc((size_t)Q + 1, sizeof(double));
	int *tau = (int *)calloc((size_t)Q + 3, sizeof(int)), *cptsout = (int *)calloc((size_t)Q + 1, sizeof(int));
	double oldmax = 1.7E+308, nullc, maxout;
	int q, p, i, j, whichout, st, end, op, stop, k, a, b;
	if (!lambda || !likeout || !tau || !cptsout) { free(lambda); free(likeout); free(tau); free(cptsout); return 909; }
	for (i = 0; i < n; i++) lambda[i] = -INFINITY;
	tau[0] = 0; tau[1] = n;
	for (q = 0; q < Q; q++) {
		for (p = 0; p < n; p++) lambda[p] = 0;
		i = 1; st = tau[0] + 1; end = tau[1];
		nullc = (-0.5) * cf(s1[end] - s1[st - 1], s2[end] - s2[st - 1], s3[end] - s3[st - 1], end - st + 1);
		for (j = 2; j < (n - 2); j++) {
			if (j == end) {
				st = end + 1; i = i + 1; end = tau[i];
				nullc = (-0.5) * cf(s1[end] - s1[st - 1], s2[end] - s2[st - 1], s3[end] - s3[st - 1], end - st + 1);
			}
			else if (((j - st) >= minseglen) && ((end - j) >= minseglen)) {
				lambda[j] = ((-0.5) * cf(s1[j] - s1[st - 1], s2[j] - s2[st - 1], s3[j] - s3[st - 1], j - st + 1)) +
					((-0.5) * cf(s1[end] - s1[j], s2[end] - s2[j], s3[end] - s3[j], end - j)) - nullc;
			}
		}
		maxout = lambda[0]; whichout = 0;
		for (i = 1; i < n; i++) if (lambda[i] > maxout) { maxout = lambda[i]; whichout = i; }
		cptsout[q] = whichout;
		likeout[q] = (oldmax <= maxout) ? oldmax : maxout;
		oldmax = likeout[q];
		tau[q + 2] = whichout;
		/* order_vec(tau, q + 3): bubble sort */
		for (a = 0; a < q + 3; a++) for (b = 1; b < (q + 3 - a); b++) if (tau[b - 1] > tau[b]) { int t = tau[b - 1]; tau[b - 1] = tau[b]; tau[b] = t; }
	}
	stop = 0; op = 0;
	while (stop == 0 && op < Q) { if ((2 * likeout[op]) >= pen) op++; else stop = 1; }
	for (k = 0; k < op; k++) cp[k] = cptsout[k];
	for (a = 1; a < op; a++) { int v = cp[a]; b = a - 1; while (b >= 0 && cp[b] > v) { cp[b + 1] = cp[b]; b--; } cp[b + 1] = v; }
	cp[op] = n;
	*ncp = op + 1;
	free(lambda); free(likeout); free(tau); free(cptsout);
	return 0;
}

/* at most one change (changepoint's single.mean.norm / single.var.norm with penalty MBIC); mean = 1, variance = 2.
   x = the data, mu = its mean (variance). Returns the change point (n = none). */
static int cpt_amoc(const double *x, int n, int moment, int minseglen, int *cp, int *ncp) {
	double *y = (double *)calloc((size_t)n + 1, sizeof(double)), *y2 = (double *)calloc((size_t)n + 1, sizeof(double));
	double nullc, tmin = 0.0, taulike, pen = 3.0 * log((double)n), alt, teststat, mu = 0.0;
	int t, found = 0, tau = 0, ntau = n - minseglen + 1 - minseglen + 1;
	if (!y || !y2) { free(y); free(y2); return 909; }
	if (moment == 2) { for (t = 0; t < n; t++) mu += x[t]; mu /= (double)n; }
	for (t = 0; t < n; t++) {
		if (moment == 1) { y[t + 1] = y[t] + x[t]; y2[t + 1] = y2[t] + x[t] * x[t]; }
		else y[t + 1] = y[t] + (x[t] - mu) * (x[t] - mu);
	}
	/* (cumsum is a running sum in R; the same order) */
	if (moment == 1) nullc = y2[n] - y[n] * y[n] / n;
	else nullc = n * log(y[n] / n);
	for (t = 0; t < ntau; t++) {
		int ts = minseglen + t;
		double v;
		if (moment == 1) v = y2[ts] - y[ts] * y[ts] / ts + (y2[n] - y2[ts]) - ((y[n] - y[ts]) * (y[n] - y[ts])) / (n - ts);
		else {
			double s1 = y[ts] / ts, sn = (y[n] - y[ts]) / (n - ts);
			if (s1 <= 0) s1 = 1e-10;
			if (sn <= 0) sn = 1e-10;
			v = ts * log(s1) + (n - ts) * log(sn);
		}
		if (isnan(v)) continue;
		if (!found || v < tmin) { tmin = v; found = 1; tau = ts; }
	}
	free(y); free(y2);
	taulike = tmin;
	alt = taulike + log((double)tau) + log((double)(n - tau + 1));
	if (!found || !isfinite(alt)) teststat = nullc; else teststat = nullc - alt;
	if (found && teststat >= pen) { cp[0] = tau; cp[1] = n; *ncp = (tau == n) ? 1 : 2; }
	else { cp[0] = n; *ncp = 1; }
	return 0;
}

/* examineOutliers: kind 1 Hampel, 2 IQR, 3 Top; returns the (truncated) indices of the outliers in out. err: 1 window not usable,
   2 parameter not usable */
static int ser_outliers(const double *x, nwint n, int kind, double par, double window, int hasw, nwint *out, nwint *nout, int *err) {
	nwint k, m = 0;
	double *buf = (double *)malloc(((size_t)n + 1) * sizeof(double));
	*err = 0;
	if (!buf) return 909;
	if (kind == 3) {
		/* order(x)[1:parameter] (stable) */
		nwint *idx = (nwint *)malloc(((size_t)n + 1) * sizeof(nwint)), cnt;
		nwint a, b;
		if (!idx) { free(buf); return 909; }
		cnt = (par < 1.0) ? 1 : (nwint)floor(par);
		if (cnt > n) { *err = 2; free(idx); free(buf); return 0; }
		for (k = 0; k < n; k++) idx[k] = k;
		/* stable insertion of the cnt smallest: partial selection sort with first-index tie rule */
		for (a = 0; a < cnt; a++) {
			nwint best = a;
			for (b = a + 1; b < n; b++) if (x[idx[b]] < x[idx[best]] || (x[idx[b]] == x[idx[best]] && idx[b] < idx[best])) best = b;
			{ nwint t = idx[a]; idx[a] = idx[best]; idx[best] = t; }
			/* keep the remaining elements in their original relative order is not needed: ties are resolved by index */
		}
		for (a = 0; a < cnt; a++) out[m++] = idx[a] + 1;
		free(idx);
	}
	else if (kind == 2) {
		double med, q1, q3, iqr, lim;
		memcpy(buf, x, (size_t)n * sizeof(double));
		med = med_sorted(buf, n);   /* buf is sorted now */
		q1 = quant7(buf, n, 0.25); q3 = quant7(buf, n, 0.75);
		iqr = q3 - q1;
		lim = med - (par / 2.0) * iqr;
		for (k = 0; k < n; k++) if (x[k] < lim) out[m++] = k + 1;
	}
	else {
		double w = hasw ? window : ((double)n / 2.0) - 1.0, from = w + 1.0, to = (double)n - w, L = 1.4826, i;
		nwint nseq, q;
		if (w < 0.0 || from > to + 1e-10) { *err = 1; free(buf); return 0; }
		nseq = (nwint)floor(to - from + 1e-10) + 1;
		for (q = 0; q < nseq; q++) {
			double lo, hi, x0, S0;
			nwint cnt = 0, t, nt;
			i = from + (double)q;
			lo = i - w; hi = i + w;
			nt = (nwint)floor(hi - lo + 1e-10) + 1;
			for (t = 0; t < nt; t++) {
				nwint ix = (nwint)(lo + (double)t);
				if (ix >= 1 && ix <= n) buf[cnt++] = x[ix - 1];
				else if (ix > n) { *err = 1; free(buf); return 0; }
			}
			x0 = med_sorted(buf, cnt);
			/* buf is sorted: recompute the absolute deviations from a fresh copy of the segment */
			for (t = 0, cnt = 0; t < nt; t++) {
				nwint ix = (nwint)(lo + (double)t);
				if (ix >= 1 && ix <= n) buf[cnt++] = fabs(x[ix - 1] - x0);
			}
			S0 = L * med_sorted(buf, cnt);
			{
				nwint xi = (nwint)i;
				if (xi >= 1 && xi <= n && fabs(x[xi - 1] - x0) > par * S0) out[m++] = xi;
			}
		}
	}
	*nout = m;
	free(buf);
	return 0;
}

/* examineChangepoints: moment 1 mean, 2 variance; method 1 PELT, 2 AMOC, 3 BinSeg; minseglen as in goldfish's window.
   The result is goldfish's cpt.pts after its de-duplication on the event times (the last of events at one time stays).
   err: 1 minseglen too large for the data, 2 too few observations, 3 Q too large */
static int ser_changepoints(const double *x, const double *tim, int n, int method, int moment, int minseglen, int *out, int *nout, int *err) {
	double *s1 = NULL, *s2 = NULL, *s3 = NULL, mu = 0.0, pen;
	int *cp = (int *)calloc((size_t)n + 3, sizeof(int)), ncp = 0, i, rc = 0, Q = 5;
	CptCost cf = (moment == 1) ? cost_mean_mbic : cost_var_mbic;
	*err = 0;
	if (!cp) return 909;
	if (moment == 2 && minseglen < 2) minseglen = 2;
	if (minseglen < 1) minseglen = 1;
	if (n < 2 || (moment == 2 && n < 4)) { *err = 2; free(cp); return 0; }
	if (n < 2 * minseglen) { *err = 1; free(cp); return 0; }
	if (method == 3 && Q > ((double)n / minseglen) + 1) { *err = 3; free(cp); return 0; }
	if (method == 3 && Q > n) { *err = 3; free(cp); return 0; }
	if (method == 2) {
		rc = cpt_amoc(x, n, moment, minseglen, cp, &ncp);
	}
	else {
		s1 = (double *)calloc((size_t)n + 1, sizeof(double)); s2 = (double *)calloc((size_t)n + 1, sizeof(double)); s3 = (double *)calloc((size_t)n + 1, sizeof(double));
		if (!s1 || !s2 || !s3) { free(s1); free(s2); free(s3); free(cp); return 909; }
		for (i = 0; i < n; i++) mu += x[i];
		mu /= (double)n;
		for (i = 0; i < n; i++) { s1[i + 1] = s1[i] + x[i]; s2[i + 1] = s2[i] + x[i] * x[i]; s3[i + 1] = s3[i] + (x[i] - mu) * (x[i] - mu); }
		pen = 3.0 * log((double)n);
		if (method == 1) rc = cpt_pelt(s1, s2, s3, n, pen, minseglen, cf, cp, &ncp);
		else rc = cpt_binseg(s1, s2, s3, n, pen, minseglen, Q, cf, cp, &ncp);
		free(s1); free(s2); free(s3);
	}
	if (rc == 0) {
		/* goldfish: if the event times of the change points repeat, keep the last of each time */
		int k = 0, j, dup = 0;
		for (i = 0; i < ncp && !dup; i++) for (j = i + 1; j < ncp; j++) if (cp[i] >= 1 && cp[j] >= 1 && tim[cp[i] - 1] == tim[cp[j] - 1]) { dup = 1; break; }
		for (i = 0; i < ncp; i++) {
			int keep = 1;
			if (dup && cp[i] >= 1) for (j = i + 1; j < ncp; j++) if (cp[j] >= 1 && tim[cp[i] - 1] == tim[cp[j] - 1]) { keep = 0; break; }
			if (keep) out[k++] = cp[i];
		}
		*nout = k;
	}
	free(cp);
	return rc;
}


/* ------------------------------------------------------------------
   Plugin entry point
   ------------------------------------------------------------------ */
/* code 40: series tools; tok = the first token after the code */
static int series_call(char *tok) {
	nwint L, nd, k, nout = 0, *oi = NULL;
	int kind, method, moment, minseglen, err = 0, rc = 0, *ci = NULL, cn = 0;
	double par, window, *x = NULL, *tim = NULL;
	L = (nwint)atof(tok);
	nd = (nwint)atof(strtok(NULL, " \t"));
	kind = (int)atof(strtok(NULL, " \t"));
	method = (int)atof(strtok(NULL, " \t"));
	moment = (int)atof(strtok(NULL, " \t"));
	minseglen = (int)atof(strtok(NULL, " \t"));
	par = atof(strtok(NULL, " \t"));
	window = atof(strtok(NULL, " \t"));
	if (nd < 1) { SF_vstore(2, 1, 0.0); SF_vstore(2, 2, 2.0); return 0; }
	x = (double *)malloc((size_t)nd * sizeof(double));
	tim = (double *)malloc((size_t)nd * sizeof(double));
	oi = (nwint *)malloc(((size_t)nd + 1) * sizeof(nwint));
	ci = (int *)malloc(((size_t)nd + 3) * sizeof(int));
	if (!x || !tim || !oi || !ci) { free(x); free(tim); free(oi); free(ci); return 909; }
	for (k = 0; k < nd; k++) { ST_double v; SF_vdata(1, (ST_int)(L + 1 + k), &v); x[k] = v; SF_vdata(1, (ST_int)(L + 1 + nd + k), &v); tim[k] = v; }
	if (kind >= 1 && kind <= 3) {
		rc = ser_outliers(x, nd, kind, par, window, window >= 0.0, oi, &nout, &err);
		if (rc == 0 && err == 0) for (k = 0; k < nout; k++) SF_vstore(1, (ST_int)(L + 1 + 2 * nd + k), (ST_double)oi[k]);
	}
	else {
		rc = ser_changepoints(x, tim, (int)nd, method, moment, minseglen, ci, &cn, &err);
		nout = cn;
		if (rc == 0 && err == 0) for (k = 0; k < nout; k++) SF_vstore(1, (ST_int)(L + 1 + 2 * nd + k), (ST_double)ci[k]);
	}
	if (rc == 0) { SF_vstore(2, 1, err ? 0.0 : (ST_double)nout); SF_vstore(2, 2, (ST_double)err); }
	free(x); free(tim); free(oi); free(ci);
	return rc;
}

STDLL stata_call(int argc, char *argv[]) {
	char *argbuf, *tok;
	int code, todo, Pt, k, rc = 0, overflow = 0, postflags = 0;
	PredOut po;
	nwint L, nrecmax = 0;
	double *theta = NULL, *g = NULL, *H = NULL, ll = 0.0;
	Spec sp;

	if (argc < 1) { SF_error("dynam_sim: missing argument string\n"); return 198; }
	argbuf = (char *)malloc(strlen(argv[0]) + 1);
	if (!argbuf) return 909;
	strcpy(argbuf, argv[0]);
	tok = strtok(argbuf, " \t");
	code = tok ? (int)atof(tok) : 0;
	if (code == 99) {
		free(argbuf);
		SF_vstore(2, 1, (ST_double)DYNAM_PROTOCOL);
		return 0;
	}
	if (code == 40) {
		char *t1 = strtok(NULL, " \t");
		int rc40 = t1 ? series_call(t1) : 198;
		free(argbuf);
		return rc40;
	}
	if (code != 20 && code != 30) {
		free(argbuf);
		SF_error("dynam_sim: unknown algorithm code (plugin and package versions do not match)\n");
		return 198;
	}
	L = (nwint)atof(strtok(NULL, " \t"));
	todo = (int)atof(strtok(NULL, " \t"));
	if (code == 30) { postflags = todo; todo = 0; }
	Pt = (int)atof(strtok(NULL, " \t"));
	theta = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	g = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	H = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * (size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	if (!theta || !g || !H) { free(argbuf); free(theta); free(g); free(H); return 909; }
	for (k = 0; k < Pt; k++) theta[k] = atof(strtok(NULL, " \t"));
	if (code == 30 && (postflags & 1)) { char *tc = strtok(NULL, " \t"); nrecmax = tc ? (nwint)atof(tc) : 0; }
	free(argbuf);

	rc = spec_read(&sp, L);
	if (rc == 0) {
		int expectP = sp.P + (sp.model == M_RATEINT ? 1 : 0);
		if (expectP != Pt) { SF_error("dynam_sim: parameter count does not match the model\n"); rc = 198; }
	}
	if (rc == 0 && code == 30) {
		memset(&po, 0, sizeof(po));
		po.write = (postflags & 1) != 0;
		po.cap = nrecmax;
		po.pstride = (sp.model == M_COORD) ? sp.n * (sp.n - 1) / 2 : sp.n;
		if (po.write) {
			nwint i2;
			po.ll = (double *)malloc((size_t)nrecmax * sizeof(double));
			po.typ = (double *)malloc((size_t)nrecmax * sizeof(double));
			po.tim = (double *)malloc((size_t)nrecmax * sizeof(double));
			if (postflags & 2) po.pr = (double *)malloc((size_t)sp.N * (size_t)po.pstride * sizeof(double));
			if (!po.ll || !po.typ || !po.tim || ((postflags & 2) && !po.pr)) rc = 909;
			(void)i2;
		}
		g_po = (rc == 0) ? &po : NULL;
	}
	if (rc == 0) {
		switch (sp.model) {
		case M_CHOICE:
		case M_RATE:
			rc = eval_softmax(&sp, theta, todo, &ll, g, H);
			break;
		case M_RATEINT:
			rc = eval_rateint(&sp, theta, todo, &ll, g, H, &overflow);
			break;
		case M_COORD:
			rc = eval_coord(&sp, theta, todo, &ll, g);
			break;
		default:
			SF_error("dynam_sim: unknown model\n");
			rc = 198;
		}
	}
	if (code == 30) {
		g_po = NULL;
		if (rc == 0 && !overflow) {
			nwint i2;
			SF_vstore(2, 2, (ST_double)po.nrec);
			SF_vstore(2, 3, (ST_double)sp.N);
			SF_vstore(2, 4, (ST_double)po.pstride);
			if (po.write) {
				for (i2 = 0; i2 < po.nrec; i2++) {
					SF_vstore(1, (ST_int)(L + 1 + i2), po.ll[i2]);
					SF_vstore(1, (ST_int)(L + 1 + po.nrec + i2), po.typ[i2]);
					SF_vstore(1, (ST_int)(L + 1 + 2 * po.nrec + i2), po.tim[i2]);
				}
				if (po.pr) {
					size_t tot = (size_t)sp.N * (size_t)po.pstride, q2;
					for (q2 = 0; q2 < tot; q2++) SF_vstore(1, (ST_int)(L + 1 + 3 * po.nrec + q2), po.pr[q2]);
				}
			}
		}
		free(po.ll); free(po.typ); free(po.tim); free(po.pr);
		if (rc == 0) SF_vstore(2, 1, overflow ? SV_missval : ll);
		free(sp.blob);
		free(theta); free(g); free(H);
		return rc;
	}
	if (rc == 0) {
		SF_vstore(2, 1, overflow ? SV_missval : ll);
		for (k = 0; k < Pt; k++) SF_vstore(2, (ST_int)(2 + k), overflow ? SV_missval : g[k]);
		if (todo >= 2 && sp.model != M_COORD) for (k = 0; k < Pt * Pt; k++) SF_vstore(2, (ST_int)(2 + Pt + k), overflow ? SV_missval : H[k]);
	}
	free(sp.blob);
	free(theta); free(g); free(H);
	return rc;
}
