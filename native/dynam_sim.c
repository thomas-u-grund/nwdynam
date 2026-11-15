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

	Blob (doubles; L = total length):
	  [0] protocol  [1] model (1 choice, 2 rate, 3 rate with intercept, 4 coordination)  [2] n  [3] N  [4] n1 (mode-1
	  actors of a two-mode network, 0 = one-mode)  [5] M (number of effect slots)  [6] opportunities given  [7] tie() given
	  [8..11]  windows inertia, recip, indeg, outdeg (>= 1e300 = none)
	  [12..15] weighted flags inertia, recip, indeg, outdeg
	  [16..35] active flags of the effect slots (choice 17, rate 5, coordination 12)
	  [36..] events: N rows of (sender, receiver, time, value, last-of-group flag), sorted by time (the caller sorts)
	  then 7 covariate rows of n values (same, diff, sim, alter/ego, tertius, egoalterint1, egoalterint2),
	  then the N x n opportunity matrix (only when [6] = 1), then the n x n matrix of tie() (only when [7] = 1).

	The statistics are those of the Mata reference evaluators (cscripts/ref_dynam_mata.do), which agree with goldfish.
	All counts are kept up to date in O(1) or O(n) per tie change instead of being recomputed at every event.
*/

#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef int64_t nwint;

#define DYNAM_PROTOCOL 2
#define NOWIN 1e300
#define HDR 36
#define M_CHOICE 1
#define M_RATE 2
#define M_RATEINT 3
#define M_COORD 4

typedef struct {
	int model;
	nwint n, N, n1;
	int M, hasopp, hastie;
	double win[4], wgt[4];
	double act[20];
	double *ev;        /* N x 5 */
	double *cov[7];    /* each n */
	double *opp;       /* N x n or NULL */
	double *tie;       /* n x n or NULL */
	double *blob;
	int eff[20], P;    /* active effect slots (1-based) and their number */
} Spec;

typedef struct {
	nwint n;
	double *tv;                        /* tie values, n x n */
	char *Ab;                          /* tv > 0 */
	double *inS, *outS, *inV, *outV;   /* column / row sums of sign(tv) and of tv */
	double *inA, *outA;                /* column / row counts of Ab */
	double *W[4];                      /* windowed tie values (slot 0 inertia, 1 recip, 2 indeg, 3 outdeg) */
	double *wInS, *wInV, *wOutS, *wOutV;
	nwint ptr[4];
	int32_t *AA, *AtA, *AAt;
	double *nt, *tnum;
	int sym;
} State;

static double sgn(double x) { return (double)((x > 0.0) - (x < 0.0)); }

static void state_free(State *st) {
	int q;
	free(st->tv); free(st->Ab); free(st->inS); free(st->outS); free(st->inV); free(st->outV);
	free(st->inA); free(st->outA);
	for (q = 0; q < 4; q++) free(st->W[q]);
	free(st->wInS); free(st->wInV); free(st->wOutS); free(st->wOutV);
	free(st->AA); free(st->AtA); free(st->AAt); free(st->nt); free(st->tnum);
	memset(st, 0, sizeof(State));
}

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
	st->tv = (double *)calloc(nn, sizeof(double));
	st->Ab = (char *)calloc(nn, 1);
	st->inS = (double *)calloc(un, sizeof(double)); st->outS = (double *)calloc(un, sizeof(double));
	st->inV = (double *)calloc(un, sizeof(double)); st->outV = (double *)calloc(un, sizeof(double));
	st->inA = (double *)calloc(un, sizeof(double)); st->outA = (double *)calloc(un, sizeof(double));
	if (!st->tv || !st->Ab || !st->inS || !st->outS || !st->inV || !st->outV || !st->inA || !st->outA) return 909;
	for (q = 0; q < 4; q++) {
		int use = (ch && q < 2) || (ra && q >= 2);
		if (use && sp->win[q] < NOWIN) {
			st->W[q] = (double *)calloc(nn, sizeof(double));
			if (!st->W[q]) return 909;
		}
	}
	if (ra && sp->win[2] < NOWIN) { st->wInS = (double *)calloc(un, sizeof(double)); st->wInV = (double *)calloc(un, sizeof(double)); if (!st->wInS || !st->wInV) return 909; }
	if (ra && sp->win[3] < NOWIN) { st->wOutS = (double *)calloc(un, sizeof(double)); st->wOutV = (double *)calloc(un, sizeof(double)); if (!st->wOutS || !st->wOutV) return 909; }
	if (ch) {
		if (has_eff(sp, 9) || has_eff(sp, 10)) st->AA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 11) || has_eff(sp, 13)) st->AtA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 12)) st->AAt = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 14)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 16)) st->tnum = (double *)calloc(un, sizeof(double));
		if (((has_eff(sp, 9) || has_eff(sp, 10)) && !st->AA) || ((has_eff(sp, 11) || has_eff(sp, 13)) && !st->AtA) ||
			(has_eff(sp, 12) && !st->AAt) || (has_eff(sp, 14) && !st->nt) || (has_eff(sp, 16) && !st->tnum)) return 909;
	}
	else if (ra) {
		if (has_eff(sp, 4)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 5)) st->tnum = (double *)calloc(un, sizeof(double));
		if ((has_eff(sp, 4) && !st->nt) || (has_eff(sp, 5) && !st->tnum)) return 909;
	}
	else if (co) {
		/* coordination: the products are recomputed per event from Ab where needed; nodeTrans and tertius are kept */
		if (has_eff(sp, 7)) st->nt = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 9)) st->tnum = (double *)calloc(un, sizeof(double));
		if (has_eff(sp, 8)) st->AA = (int32_t *)calloc(nn, sizeof(int32_t));
		if (has_eff(sp, 10)) st->AtA = (int32_t *)calloc(nn, sizeof(int32_t));
		if ((has_eff(sp, 7) && !st->nt) || (has_eff(sp, 9) && !st->tnum) || (has_eff(sp, 8) && !st->AA) || (has_eff(sp, 10) && !st->AtA)) return 909;
	}
	return 0;
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

static void wupd(State *st, int q, nwint a, nwint b, double dv) {
	nwint n = st->n;
	double old = st->W[q][a * n + b], nw = old + dv;
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
		if (st->tnum) st->tnum[b] += (na - oa) * sp->cov[4][a];
		st->Ab[a * n + b] = (char)na;
	}
	st->tv[a * n + b] = nw;
	for (q = 0; q < 4; q++) if (st->W[q]) wupd(st, q, a, b, v);
}

static void apply_event(const Spec *sp, State *st, nwint k) {
	nwint a = (nwint)sp->ev[5 * k] - 1, b = (nwint)sp->ev[5 * k + 1] - 1;
	double v = sp->ev[5 * k + 3];
	apply_cell(sp, st, a, b, v);
	if (st->sym) apply_cell(sp, st, b, a, v);
}

/* windows: drop the events of earlier groups (rows < gs) that are older than the window at time tcur */
static void expire(const Spec *sp, State *st, nwint gs, double tcur) {
	int q;
	nwint n = st->n;
	for (q = 0; q < 4; q++) {
		if (!st->W[q]) continue;
		while (st->ptr[q] < gs) {
			nwint p = st->ptr[q];
			if (tcur - sp->ev[5 * p + 2] <= sp->win[q]) break;
			wupd(st, q, (nwint)sp->ev[5 * p] - 1, (nwint)sp->ev[5 * p + 1] - 1, -sp->ev[5 * p + 3]);
			st->ptr[q]++;
		}
	}
	(void)n;
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
		case 4: for (j = 0; j < n; j++) x[j] = (sp->cov[0][j] == sp->cov[0][s]); break;
		case 5: for (j = 0; j < n; j++) x[j] = fabs(sp->cov[1][j] - sp->cov[1][s]); break;
		case 6: for (j = 0; j < n; j++) x[j] = -fabs(sp->cov[2][j] - sp->cov[2][s]); break;
		case 7: for (j = 0; j < n; j++) x[j] = sp->cov[3][j]; break;
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
		case 15: for (j = 0; j < n; j++) x[j] = sp->cov[5][s] * sp->cov[6][j]; break;
		case 16: for (j = 0; j < n; j++) x[j] = (st->inA[j] > 0.0) ? st->tnum[j] / st->inA[j] : 0.0; break;
		case 17: for (j = 0; j < n; j++) x[j] = sp->tie[s * n + j]; break;
		}
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
		case 3: for (j = 0; j < n; j++) x[j] = sp->cov[3][j]; break;
		case 4: for (j = 0; j < n; j++) x[j] = st->nt[j]; break;
		case 5: for (j = 0; j < n; j++) x[j] = (st->inA[j] > 0.0) ? st->tnum[j] / st->inA[j] : 0.0; break;
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
			if (!choice) rate_stats(sp, &st, X, 0);
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
		mx = lp[0];
		for (j = 1; j < n; j++) if (lp[j] > mx) mx = lp[j];
		sumP = 0.0;
		for (j = 0; j < n; j++) { Pr[j] = exp(lp[j] - mx); sumP += Pr[j]; }
		lrsum = mx + log(sumP);
		L += lp[choice ? r : s] - lrsum;
		for (j = 0; j < n; j++) Pr[j] /= sumP;
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
static int eval_rateint(const Spec *sp, const double *theta, int todo, double *ll, double *g, double *H, int *overflow) {
	State st;
	nwint n = sp->n, N = sp->N, i, j, gs = 0, nrisk = (sp->n1 > 0 ? sp->n1 : sp->n);
	int Pt = sp->P + 1, k, l, rc = 0;
	double *X = NULL, *haz = NULL, *Sx = NULL, *Sxx = NULL, L = 0.0, tprev = 0.0, w, S;
	int fresh = 1;

	*overflow = 0;
	if ((rc = state_init(sp, &st)) != 0) { state_free(&st); return rc; }
	X = (double *)malloc((size_t)Pt * (size_t)n * sizeof(double));
	haz = (double *)malloc((size_t)n * sizeof(double));
	Sx = (double *)malloc((size_t)Pt * sizeof(double));
	Sxx = (double *)malloc((size_t)Pt * (size_t)Pt * sizeof(double));
	if (!X || !haz || !Sx || !Sxx) { rc = 909; goto done; }
	for (k = 0; k < Pt; k++) g[k] = 0.0;
	if (todo >= 2) for (k = 0; k < Pt * Pt; k++) H[k] = 0.0;
	/* the hazard sums (S, Sx, Sxx) depend on the state only: computed once per group, used by its first event (w > 0) */
	for (i = 0; i < N; i++) {
		nwint s = (nwint)sp->ev[5 * i] - 1;
		double tcur = sp->ev[5 * i + 2], lps = 0.0;
		if (fresh) {
			expire(sp, &st, gs, tcur);
			rate_stats(sp, &st, X, 1);
		}
		w = (i > 0) ? tcur - tprev : 0.0;
		for (k = 0; k < Pt; k++) lps += theta[k] * X[(size_t)k * (size_t)n + s];
		if (lps > 700.0) { *overflow = 1; goto done; }
		L += lps;
		for (k = 0; k < Pt; k++) g[k] += X[(size_t)k * (size_t)n + s];
		if (w != 0.0) {
			S = 0.0;
			for (k = 0; k < Pt; k++) Sx[k] = 0.0;
			if (todo >= 2) for (k = 0; k < Pt * Pt; k++) Sxx[k] = 0.0;
			for (j = 0; j < nrisk; j++) {
				double lp = 0.0, ej;
				for (k = 0; k < Pt; k++) lp += theta[k] * X[(size_t)k * (size_t)n + j];
				if (lp > 700.0) { *overflow = 1; goto done; }
				ej = exp(lp);
				haz[j] = ej;
				S += ej;
				for (k = 0; k < Pt; k++) {
					double xk = ej * X[(size_t)k * (size_t)n + j];
					Sx[k] += xk;
					if (todo >= 2) for (l = k; l < Pt; l++) Sxx[k * Pt + l] += xk * X[(size_t)l * (size_t)n + j];
				}
			}
			L -= w * S;
			for (k = 0; k < Pt; k++) g[k] -= w * Sx[k];
			if (todo >= 2) for (k = 0; k < Pt; k++) for (l = k; l < Pt; l++) {
				double h = -w * Sxx[k * Pt + l];
				H[k * Pt + l] += h;
				if (l != k) H[l * Pt + k] += h;
			}
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
	*ll = L;
done:
	free(X); free(haz); free(Sx); free(Sxx);
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
static double coord_S(const Spec *sp, const State *st, int e, nwint a, nwint b, const int32_t *AAtA) {
	nwint n = sp->n;
	switch (e) {
	case 1: return (double)st->Ab[a * n + b];
	case 2: return st->inS[b];
	case 3: return (double)(sp->cov[0][a] == sp->cov[0][b]);
	case 4: return fabs(sp->cov[1][a] - sp->cov[1][b]);
	case 5: return -fabs(sp->cov[2][a] - sp->cov[2][b]);
	case 6: return sp->cov[3][b];
	case 7: return st->nt[b];
	case 8: return (double)st->AA[a * n + b];
	case 9: return (st->inA[b] > 0.0) ? st->tnum[b] / st->inA[b] : 0.0;
	case 10: {
		double A = (double)st->Ab[a * n + b];
		return (double)AAtA[a * n + b] - A * (st->outA[a] + (double)st->AtA[b * n + b] - A);
	}
	case 11: return sp->cov[5][a] * sp->cov[6][b];
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
		for (a = 0; a < n; a++) for (b = 0; b < n; b++) {
			double v = 0.0;
			for (q = 0; q < P; q++) v += theta[q] * coord_S(sp, &st, sp->eff[q], a, b, AAtA);
			stm[a * n + b] = v;
		}
		for (a = 0; a < n; a++) {
			double rmax = stm[a * n], rs = 0.0;
			for (b = 1; b < n; b++) if (stm[a * n + b] > rmax) rmax = stm[a * n + b];
			for (b = 0; b < n; b++) { double e = (a == b) ? 0.0 : exp(stm[a * n + b] - rmax); p[a * n + b] = e; rs += e; }
			for (b = 0; b < n; b++) p[a * n + b] /= rs;
		}
		tot = 0.0;
		for (a = 0; a < n; a++) for (b = 0; b < n; b++) { double v = p[a * n + b] * p[b * n + a]; Pn[a * n + b] = v; tot += v; }
		Z = tot / 2.0;
		for (a = 0; a < n; a++) { double rs = 0.0; for (b = 0; b < n; b++) { Pn[a * n + b] /= Z; rs += Pn[a * n + b]; } rowPn[a] = rs; }
		L += log(Pn[s * n + r]);
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
	for (k = 0; k < 20; k++) sp->act[k] = sp->blob[16 + k];
	sp->P = 0;
	for (k = 0; k < sp->M && k < 20; k++) if (sp->act[k] != 0.0) sp->eff[sp->P++] = k + 1;
	pos = HDR;
	sp->ev = sp->blob + pos; pos += 5 * sp->N;
	for (k = 0; k < 7; k++) { sp->cov[k] = sp->blob + pos; pos += sp->n; }
	if (sp->hasopp) { sp->opp = sp->blob + pos; pos += sp->N * sp->n; }
	if (sp->hastie) { sp->tie = sp->blob + pos; pos += sp->n * sp->n; }
	if (pos != L) return 198;
	return 0;
}

/* ------------------------------------------------------------------
   Plugin entry point
   ------------------------------------------------------------------ */
STDLL stata_call(int argc, char *argv[]) {
	char *argbuf, *tok;
	int code, todo, Pt, k, rc = 0, overflow = 0;
	nwint L;
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
	if (code != 20) {
		free(argbuf);
		SF_error("dynam_sim: unknown algorithm code (plugin and package versions do not match)\n");
		return 198;
	}
	L = (nwint)atof(strtok(NULL, " \t"));
	todo = (int)atof(strtok(NULL, " \t"));
	Pt = (int)atof(strtok(NULL, " \t"));
	theta = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	g = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	H = (double *)malloc((size_t)(Pt > 0 ? Pt : 1) * (size_t)(Pt > 0 ? Pt : 1) * sizeof(double));
	if (!theta || !g || !H) { free(argbuf); free(theta); free(g); free(H); return 909; }
	for (k = 0; k < Pt; k++) theta[k] = atof(strtok(NULL, " \t"));
	free(argbuf);

	rc = spec_read(&sp, L);
	if (rc == 0) {
		int expectP = sp.P + (sp.model == M_RATEINT ? 1 : 0);
		if (expectP != Pt) { SF_error("dynam_sim: parameter count does not match the model\n"); rc = 198; }
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
	if (rc == 0) {
		SF_vstore(2, 1, overflow ? SV_missval : ll);
		for (k = 0; k < Pt; k++) SF_vstore(2, (ST_int)(2 + k), overflow ? SV_missval : g[k]);
		if (todo >= 2 && sp.model != M_COORD) for (k = 0; k < Pt * Pt; k++) SF_vstore(2, (ST_int)(2 + Pt + k), overflow ? SV_missval : H[k]);
	}
	free(sp.blob);
	free(theta); free(g); free(H);
	return rc;
}
