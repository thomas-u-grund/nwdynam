* dev/dynam_unit2_rate_crosscheck.do -- fits nwdynam's own Unit 2
* rate-submodel engine (unw_dynam.do's DynamFitRateUnit1()) on the EXACT
* SAME event data (dev/dynam_crosscheck_events.csv, the real goldfish
* Social_Evolution `calls` dataset) the real
* goldfish::estimate(model="DyNAM", subModel="rate") run fit - matching
* dev/dynam_unit1_choice_crosscheck.do's own role for the choice
* sub-model and nwrem's own dev/rem_relevent_crosscheck.do precedent.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439 events,
* `callsDependent ~ indeg(callNetwork) + outdeg(callNetwork)`,
* model="DyNAM", subModel="rate", NO intercept - see
* dev/dynam_unit1_crosscheck.R's own real console output, captured
* 2026-09-02):
*   indeg  =  0.551445  (SE 0.066344)
*   outdeg =  0.263784  (SE 0.028386)
*   logLik = -1750.9

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
bname = "dynam_rate_b"
vname = "dynam_rate_V"
llname = "dynam_rate_ll"
DynamFitRateUnit1(eventmat, n_actors, bname, vname, llname)
end

di as text "{hline 60}"
di as text "nwdynam rate sub-model (Unit 2, no intercept) fit on the SAME event data as"
di as text "the real goldfish::estimate(model=" `"""' "DyNAM" `"""' ", subModel=" `"""' "rate" `"""' ") run above:"
di as text "  indeg  = " %9.6f el(dynam_rate_b,1,1) "   (goldfish:  0.551445)"
di as text "  outdeg = " %9.6f el(dynam_rate_b,1,2) "   (goldfish:  0.263784)"
di as text "  logLik = " %9.3f dynam_rate_ll        "   (goldfish: -1750.9)"

scalar reldiff_indeg  = abs(el(dynam_rate_b,1,1) - 0.551445)
scalar reldiff_outdeg = abs(el(dynam_rate_b,1,2) - 0.263784)
scalar reldiff_ll     = abs(dynam_rate_ll - (-1750.9))
di as text "{hline 60}"
di as text "abs difference vs. real goldfish: indeg=" %9.6f reldiff_indeg ///
	"  outdeg=" %9.6f reldiff_outdeg "  logLik=" %9.4f reldiff_ll

* Same bar as dev/dynam_unit1_choice_crosscheck.do / dev/rem_relevent_crosscheck.do:
* both implementations climb the identical likelihood surface (see
* dev/dynam_unit2_rate_toy_crosscheck.R for the derivation that
* established this is the right surface) and should land at the
* identical optimum.
assert reldiff_indeg < 0.001
assert reldiff_outdeg < 0.001
assert reldiff_ll < 0.05

di as result "MATCH CONFIRMED (agreement to better than 1e-3 on both coefficients, 5e-2 on logLik)"
