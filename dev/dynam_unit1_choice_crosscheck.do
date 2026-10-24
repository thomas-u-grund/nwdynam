* dev/dynam_unit1_choice_crosscheck.do -- fits nwdynam's own v1 Unit 1
* choice-submodel engine (unw_dynam.do's DynamFitUnit1()) on the EXACT
* SAME event data (dev/dynam_crosscheck_events.csv, the real goldfish
* Social_Evolution `calls` dataset, exported by dev/dynam_unit1_crosscheck.R)
* the real goldfish::estimate(model="DyNAM", subModel="choice") run fit -
* matching nwrem's own development discipline
* (dev/rem_relevent_crosscheck.do) of checking a from-scratch engine
* against a real external reference implementation on IDENTICAL data,
* not a self-consistency check alone.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439 events,
* `callsDependent ~ inertia + recip + indeg`, model="DyNAM",
* subModel="choice" - see dev/dynam_unit1_crosscheck.R's own real console
* output, captured 2026-09-02):
*   inertia =  5.36556  (SE 0.19375)
*   recip   =  1.38944  (SE 0.16729)
*   indeg   = -0.17353  (SE 0.09097)
*   logLik  = -695.46

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
bname = "dynam_b"
vname = "dynam_V"
llname = "dynam_ll"
DynamFitUnit1(eventmat, n_actors, bname, vname, llname)
end

di as text "{hline 60}"
di as text "nwdynam choice sub-model fit on the SAME event data as the real"
di as text "goldfish::estimate(model=" `"""' "DyNAM" `"""' ", subModel=" `"""' "choice" `"""' ") run above:"
di as text "  inertia = " %9.5f el(dynam_b,1,1) "   (goldfish:  5.36556)"
di as text "  recip   = " %9.5f el(dynam_b,1,2) "   (goldfish:  1.38944)"
di as text "  indeg   = " %9.5f el(dynam_b,1,3) "   (goldfish: -0.17353)"
di as text "  logLik  = " %9.2f dynam_ll        "   (goldfish: -695.46)"

scalar reldiff_inertia = abs(el(dynam_b,1,1) - 5.36556)
scalar reldiff_recip   = abs(el(dynam_b,1,2) - 1.38944)
scalar reldiff_indeg   = abs(el(dynam_b,1,3) - (-0.17353))
scalar reldiff_ll      = abs(dynam_ll - (-695.46))
di as text "{hline 60}"
di as text "abs difference vs. real goldfish: inertia=" %9.5f reldiff_inertia ///
	"  recip=" %9.5f reldiff_recip "  indeg=" %9.5f reldiff_indeg ///
	"  logLik=" %9.4f reldiff_ll

* Both are MLE for what is claimed to be the identical log-likelihood on
* the identical data - if the effect definitions read from goldfish's own
* doc/goldfishEffects.Rmd (inertia/recip/indeg-alter, all unweighted, see
* unw_dynam.do's own header comment) were transcribed correctly, both
* implementations are climbing the identical likelihood surface and
* should land at (numerically) the identical optimum, same bar
* dev/rem_relevent_crosscheck.do already applies to nwrem vs. relevent.
assert reldiff_inertia < 0.01
assert reldiff_recip < 0.01
assert reldiff_indeg < 0.01
assert reldiff_ll < 0.05

di as result "MATCH CONFIRMED (agreement to better than 1e-2 on all three coefficients, 5e-2 on logLik)"
