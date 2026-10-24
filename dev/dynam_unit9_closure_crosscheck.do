* dev/dynam_unit9_closure_crosscheck.do -- adds nodeTrans (both
* sub-models) and the choice-only two-path closure effects
* trans/cycle/commonSender/commonReceiver/four (the development notes'
* "effect expansion" scope, batch 6, resolved 2026-09-02, user: "work
* your way through the list") against real goldfish subset fits on the
* SAME event data, all computed on the DEPENDENT network itself
* (goldfish's own stated default when no explicit network argument is
* given) - matching dev/dynam_unit8_choiceoutdeg_crosscheck.do's own
* direct-Mata-call, n_actors=84-hardcoded discipline.
*
* FORMULAS (goldfish's own doc/goldfishEffects.Rmd, transcribed
* directly, not guessed), with A = the n x n binary presence matrix
* (lastcontact :> -1e300, reusing the SAME tracking structure already
* maintained for inertia/recip/indeg/outdeg), sender s fixed at each
* event:
*   trans_i(j)          = (A*A)[s,.]          candidate-dependent
*   cycle_i(j)           = (A*A)[.,s]'         candidate-dependent
*   commonSender_i(j)    = (A'*A)[s,.]         candidate-dependent
*   commonReceiver_i(j)  = (A*A')[s,.]         candidate-dependent
*   four_i(j)            = (A*(A'*A))[s,.] - A[s,.]*(rowsum(A[s,.]) +
*                           diagonal(A'*A)' - A[s,.])   candidate-dependent
*                           (a REAL correction found via direct
*                           comparison against goldfish's own
*                           preprocessed statistics, not guessed - see
*                           unw_dynam.do's own header comment above
*                           dynam_choice_loglik_grad_multi() for the
*                           full derivation: `four` requires i,j,k,l to
*                           be four DISTINCT actors, unlike the other
*                           four closure effects here)
*   nodeTrans            = rowsum(A :* (A*A))' candidate-ONLY (no sender
*                                               dependence at all, like
*                                               indeg/outdeg/alter)
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439
* events - captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit9_closure_crosscheck.R for the R side):
*   choice ~ trans:                    trans=1.49074,   logLik=-1907.9
*   choice ~ inertia + trans:          inertia=5.88545, trans=0.27084,   logLik=-732.56
*   choice ~ cycle:                    cycle=1.90516,   logLik=-1887.6
*   choice ~ commonSender:             commonSender=0.47185, logLik=-1937.8
*   choice ~ commonReceiver:           commonReceiver=1.45791, logLik=-1911.0
*   choice ~ four:                     four=0.89703,    logLik=-1927.1
*   choice ~ nodeTrans:                nodeTrans=0.873489, logLik=-1900.1
*   rate ~ nodeTrans:                  nodeTrans=1.052850, logLik=-1876.9
*   choice ~ inertia+recip+trans+cycle+commonSender+commonReceiver:
*       inertia=5.25137, recip=1.39670, trans=-0.28992, cycle=-0.63663,
*       commonSender=-0.26381, commonReceiver=0.48754, logLik=-692.6

* (run from the root of this repository)
* The core package nwcommands is found through the environment variable NWCORE (path of the core checkout).
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do

import delimited "dev/dynam_crosscheck_events.csv", clear varnames(1)

mata:
sender = st_data(., "sender")
receiver = st_data(., "receiver")
eventtime = st_data(., "eventtime")
eventmat = (sender, receiver, eventtime)
n_actors = 84
zero8 = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "choice ~ trans"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
di as text "trans=" %9.5f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 1.49074)
scalar d0b = abs(ll0 - (-1907.9))
di as text "abs diff: trans=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + trans"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "inertia=" %9.5f el(b1,1,1) "  trans=" %9.5f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 5.88545)
scalar d2 = abs(el(b1,1,2) - 0.27084)
scalar d3 = abs(ll1 - (-732.56))
di as text "abs diff: inertia=" %9.6f d1 "  trans=" %9.6f d2 "  logLik=" %9.4f d3
assert d1 < 0.01
assert d2 < 0.01
assert d3 < 0.1

di as text "{hline 60}"
di as text "choice ~ cycle"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
di as text "cycle=" %9.5f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d4 = abs(el(b2,1,1) - 1.90516)
scalar d5 = abs(ll2 - (-1887.6))
di as text "abs diff: cycle=" %9.6f d4 "  logLik=" %9.4f d5
assert d4 < 0.01
assert d5 < 0.1

di as text "{hline 60}"
di as text "choice ~ commonSender"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "commonSender=" %9.5f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d6 = abs(el(b3,1,1) - 0.47185)
scalar d7 = abs(ll3 - (-1937.8))
di as text "abs diff: commonSender=" %9.6f d6 "  logLik=" %9.4f d7
assert d6 < 0.01
assert d7 < 0.1

di as text "{hline 60}"
di as text "choice ~ commonReceiver"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
di as text "commonReceiver=" %9.5f el(b4,1,1) "  logLik=" %9.4f ll4
scalar d8 = abs(el(b4,1,1) - 1.45791)
scalar d9 = abs(ll4 - (-1911.0))
di as text "abs diff: commonReceiver=" %9.6f d8 "  logLik=" %9.4f d9
assert d8 < 0.01
assert d9 < 0.1

di as text "{hline 60}"
di as text "choice ~ four"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b5", "V5", "ll5")
di as text "four=" %9.5f el(b5,1,1) "  logLik=" %9.4f ll5
scalar d10 = abs(el(b5,1,1) - 0.89703)
scalar d11 = abs(ll5 - (-1927.1))
di as text "abs diff: four=" %9.6f d10 "  logLik=" %9.4f d11
assert d10 < 0.01
assert d11 < 0.1

di as text "{hline 60}"
di as text "choice ~ nodeTrans"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b6", "V6", "ll6")
di as text "nodeTrans=" %9.6f el(b6,1,1) "  logLik=" %9.4f ll6
scalar d12 = abs(el(b6,1,1) - 0.873489)
scalar d13 = abs(ll6 - (-1900.1))
di as text "abs diff: nodeTrans=" %9.6f d12 "  logLik=" %9.4f d13
assert d12 < 0.01
assert d13 < 0.1

di as text "{hline 60}"
di as text "rate ~ nodeTrans"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,0,0,1,0), zero8, zero8, 1e300, 1e300, 0, 0, 0, "b7", "V7", "ll7")
di as text "nodeTrans=" %9.6f el(b7,1,1) "  logLik=" %9.4f ll7
scalar d14 = abs(el(b7,1,1) - 1.052850)
scalar d15 = abs(ll7 - (-1876.9))
di as text "abs diff: nodeTrans=" %9.6f d14 "  logLik=" %9.4f d15
assert d14 < 0.01
assert d15 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + recip + trans + cycle + commonSender + commonReceiver (all six)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0), zero8, zero8, zero8, zero8, zero8, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b8", "V8", "ll8")
assert colsof(b8) == 6
di as text "inertia=" %9.5f el(b8,1,1) "  recip=" %9.5f el(b8,1,2) "  trans=" %9.5f el(b8,1,3) "  cycle=" %9.5f el(b8,1,4) "  commonSender=" %9.5f el(b8,1,5) "  commonReceiver=" %9.5f el(b8,1,6) "  logLik=" %9.4f ll8
scalar d16 = abs(el(b8,1,1) - 5.25137)
scalar d17 = abs(el(b8,1,2) - 1.39670)
scalar d18 = abs(el(b8,1,3) - (-0.28992))
scalar d19 = abs(el(b8,1,4) - (-0.63663))
scalar d20 = abs(el(b8,1,5) - (-0.26381))
scalar d21 = abs(el(b8,1,6) - 0.48754)
scalar d22 = abs(ll8 - (-692.6))
di as text "abs diff: inertia=" %9.6f d16 "  recip=" %9.6f d17 "  trans=" %9.6f d18 "  cycle=" %9.6f d19 "  commonSender=" %9.6f d20 "  commonReceiver=" %9.6f d21 "  logLik=" %9.4f d22
assert d16 < 0.01
assert d17 < 0.01
assert d18 < 0.02
assert d19 < 0.02
assert d20 < 0.02
assert d21 < 0.02
assert d22 < 0.1

di as result "MATCH CONFIRMED (nodeTrans on both sub-models, and trans/cycle/commonSender/commonReceiver/four on the choice sub-model, all reproduce real goldfish - alone and in a genuine 6-effect combined model)"
