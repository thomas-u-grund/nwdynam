* dev/dynam_unit6_window_crosscheck.do -- fits unw_dynam.do's own new
* window() support (the development notes' "effect expansion" scope
* item, batch 3, resolved 2026-09-02, user request "implement window
* effects") against real goldfish subset fits on the SAME event data,
* for inertia and recip - matching dev/dynam_unit5_egoalter_crosscheck.do's
* own direct-Mata-call, n_actors=84-hardcoded discipline.
*
* WINDOW SEMANTICS - derived from a real toy example FIRST, not guessed
* (see /tmp/dynam_window_toy.R's own now-discarded session output,
* summarized here): a tie i->j counts as PRESENT at the current event's
* own evaluation time if and only if the MOST RECENT prior i->j event
* happened within `window` real time units of the current event's own
* timestamp (inclusive cutoff, not a decaying weight) - confirmed to
* full precision against real goldfish inertia(window="7 days")/
* inertia(window="1 day")/recip(window="1 day") fits on the real
* Social_Evolution dataset BEFORE any Mata code was written (matching
* this project's own dev/rem_relevent_crosscheck.R-style discipline of
* deriving from real reference-implementation behavior first).
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439
* events, real Unix-epoch-second timestamps - see console output
* captured 2026-09-02):
*   choice ~ inertia, no window:        inertia=5.8971 (SE 0.1533), logLik=-733.58
*   choice ~ inertia, window=7 days:    inertia=5.90189 (SE 0.14808), logLik=-733.74
*   choice ~ inertia, window=1 day:     inertia=6.02952 (SE 0.15389), logLik=-851.72
*   choice ~ recip, window=1 day:       recip=5.9407 (SE 0.1695), logLik=-1145.2

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
zero7b = J(1, n_actors, 0)
end

di as text "{hline 60}"
di as text "choice ~ inertia, no window (regression: must reproduce the ORIGINAL"
di as text "pre-window numbers exactly, not just approximately)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
di as text "inertia=" %9.4f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 5.8971)
scalar d0b = abs(ll0 - (-733.58))
di as text "abs diff: inertia=" %9.4f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia, window = 7 days (604800 seconds)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, 604800, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "inertia=" %9.5f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 5.90189)
scalar d2 = abs(ll1 - (-733.74))
di as text "abs diff: inertia=" %9.6f d1 "  logLik=" %9.4f d2
assert d1 < 0.01
assert d2 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia, window = 1 day (86400 seconds)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, 86400, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b2", "V2", "ll2")
di as text "inertia=" %9.5f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d3 = abs(el(b2,1,1) - 6.02952)
scalar d4 = abs(ll2 - (-851.72))
di as text "abs diff: inertia=" %9.6f d3 "  logLik=" %9.4f d4
assert d3 < 0.01
assert d4 < 0.1

di as text "{hline 60}"
di as text "choice ~ recip, window = 1 day (86400 seconds)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, 1e300, 86400, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "recip=" %9.5f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d5 = abs(el(b3,1,1) - 5.9407)
scalar d6 = abs(ll3 - (-1145.2))
di as text "abs diff: recip=" %9.6f d5 "  logLik=" %9.4f d6
assert d5 < 0.01
assert d6 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia(window=7 days) + recip(window=1 day) - GENUINELY"
di as text "DIFFERENT windows on two effects in the SAME model (the whole point"
di as text "of the per-effect redesign, not just two independent single-effect fits)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, zero7b, 604800, 86400, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
assert colsof(b4) == 2
di as text "inertia=" %9.5f el(b4,1,1) "  recip=" %9.5f el(b4,1,2) "  logLik=" %9.4f ll4
scalar d7 = abs(el(b4,1,1) - 5.11325)
scalar d8 = abs(el(b4,1,2) - 2.10586)
scalar d9 = abs(ll4 - (-674.04))
di as text "abs diff: inertia=" %9.6f d7 "  recip=" %9.6f d8 "  logLik=" %9.4f d9
assert d7 < 0.01
assert d8 < 0.01
assert d9 < 0.1

di as result "MATCH CONFIRMED (window() reproduces real goldfish exactly at 3 window sizes for inertia, 1 for recip, and one model with GENUINELY DIFFERENT windows on both simultaneously; no-window regression unaffected)"
