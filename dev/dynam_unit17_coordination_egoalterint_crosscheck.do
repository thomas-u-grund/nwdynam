* dev/dynam_unit17_coordination_egoalterint_crosscheck.do --
* choice_coordination effect expansion, batch 14, resolved 2026-09-02 -
* adds egoAlterInt(), verified against real goldfish subset fits on
* the SAME 6-actor/12-event toy network dev/dynam_unit14_coordination_
* crosscheck.do already uses.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit17_coordination_egoalterint_crosscheck.R for the R
* side):
*   egoAlterInt(floor,floor2) alone:      egoAlterInt=-0.03524755, logLik=-31.98356
*   inertia + egoAlterInt(floor,floor2):  inertia=-0.39844340, egoAlterInt=-0.03433753,
*                                          logLik=-31.60121

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

import delimited "dev/dynam_crosscheck_coord3_events.csv", clear varnames(1)
gen senderidx = real(substr(sender, 2, .))
gen receiveridx = real(substr(receiver, 2, .))
mata:
sender = st_data(., "senderidx")
receiver = st_data(., "receiveridx")
eventtime = st_data(., "time")
eventmat3 = (sender, receiver, eventtime)
n3 = 6
zero3 = J(1, n3, 0)
floorvec3 = (1,1,2,2,3,3)
floor2vec3 = (5,3,8,1,9,2)
end

di as text "{hline 60}"
di as text "choice_coordination ~ egoAlterInt(floor, floor2)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat3, n3, (0,0,0,0,0,0,0,0,0,0,1,0), zero3, zero3, zero3, zero3, zero3, floorvec3, floor2vec3, J(0,0,0), "b0", "V0", "ll0")
di as text "egoAlterInt=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - (-0.03524755))
scalar d0b = abs(ll0 - (-31.98356))
di as text "abs diff: egoAlterInt=" %9.6f d0a "  logLik=" %9.6f d0b
assert d0a < 0.001
assert d0b < 0.001

di as text "{hline 60}"
di as text "choice_coordination ~ inertia + egoAlterInt(floor, floor2)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat3, n3, (1,0,0,0,0,0,0,0,0,0,1,0), zero3, zero3, zero3, zero3, zero3, floorvec3, floor2vec3, J(0,0,0), "b1", "V1", "ll1")
assert colsof(b1) == 2
di as text "inertia=" %9.6f el(b1,1,1) "  egoAlterInt=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - (-0.39844340))
scalar d1b = abs(el(b1,1,2) - (-0.03433753))
scalar d1c = abs(ll1 - (-31.60121))
di as text "abs diff: inertia=" %9.6f d1a "  egoAlterInt=" %9.6f d1b "  logLik=" %9.6f d1c
assert d1a < 0.001
assert d1b < 0.001
assert d1c < 0.001

di as result "MATCH CONFIRMED (choice_coordination's egoAlterInt reproduces real goldfish exactly, alone and combined with inertia)"
