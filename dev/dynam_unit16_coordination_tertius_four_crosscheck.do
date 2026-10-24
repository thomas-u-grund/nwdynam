* dev/dynam_unit16_coordination_tertius_four_crosscheck.do --
* choice_coordination effect expansion, batch 13, resolved 2026-09-02 -
* adds tertius() and four, verified against real goldfish subset fits
* on the SAME 6-actor/12-event toy network dev/dynam_unit14_
* coordination_crosscheck.do already uses.
*
* Real goldfish reference (captured 2026-09-02, goldfish 1.6.12, see
* dev/dynam_unit16_coordination_tertius_four_crosscheck.R for the R
* side):
*   tertius(floor) alone:            tertius=0.2466349, logLik=-31.90115
*   four alone:                      four=0.2065944, logLik=-32.28807
*   inertia + tertius(floor) + four: inertia=-0.57080704, tertius=0.32457616,
*                                     four=0.01745949, logLik=-31.10879

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
end

di as text "{hline 60}"
di as text "choice_coordination ~ tertius(floor)"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat3, n3, (0,0,0,0,0,0,0,0,1,0,0,0), zero3, zero3, zero3, zero3, floorvec3, zero3, zero3, J(0,0,0), "b0", "V0", "ll0")
di as text "tertius=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.2466349)
scalar d0b = abs(ll0 - (-31.90115))
di as text "abs diff: tertius=" %9.6f d0a "  logLik=" %9.6f d0b
assert d0a < 0.001
assert d0b < 0.001

di as text "{hline 60}"
di as text "choice_coordination ~ four"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat3, n3, (0,0,0,0,0,0,0,0,0,1,0,0), zero3, zero3, zero3, zero3, zero3, zero3, zero3, J(0,0,0), "b1", "V1", "ll1")
di as text "four=" %9.6f el(b1,1,1) "  logLik=" %9.4f ll1
scalar d1a = abs(el(b1,1,1) - 0.2065944)
scalar d1b = abs(ll1 - (-32.28807))
di as text "abs diff: four=" %9.6f d1a "  logLik=" %9.6f d1b
assert d1a < 0.001
assert d1b < 0.001

di as text "{hline 60}"
di as text "choice_coordination ~ inertia + tertius(floor) + four"
di as text "{hline 60}"
mata: DynamCoordFitMulti(eventmat3, n3, (1,0,0,0,0,0,0,0,1,1,0,0), zero3, zero3, zero3, zero3, floorvec3, zero3, zero3, J(0,0,0), "b2", "V2", "ll2")
assert colsof(b2) == 3
di as text "inertia=" %9.6f el(b2,1,1) "  tertius=" %9.6f el(b2,1,2) "  four=" %9.6f el(b2,1,3) "  logLik=" %9.4f ll2
scalar d2a = abs(el(b2,1,1) - (-0.57080704))
scalar d2b = abs(el(b2,1,2) - 0.32457616)
scalar d2c = abs(el(b2,1,3) - 0.01745949)
scalar d2d = abs(ll2 - (-31.10879))
di as text "abs diff: inertia=" %9.6f d2a "  tertius=" %9.6f d2b "  four=" %9.6f d2c "  logLik=" %9.6f d2d
assert d2a < 0.001
assert d2b < 0.001
assert d2c < 0.001
assert d2d < 0.001

di as result "MATCH CONFIRMED (choice_coordination's tertius and four reproduce real goldfish exactly, alone and combined with inertia)"
