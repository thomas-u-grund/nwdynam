* dev/dynam_unit10_tertius_egoalterint_crosscheck.do -- adds tertius
* (both sub-models) and egoAlterInt (choice-only) (the development notes'
* "effect expansion" scope, batch 7, resolved 2026-09-02, continuing
* "work your way through the list") against real goldfish subset fits
* on the SAME event data - matching dev/dynam_unit9_closure_crosscheck.do's
* own direct-Mata-call, n_actors=84-hardcoded discipline.
*
* FORMULAS (goldfish's own doc/goldfishEffects.Rmd, transcribed
* directly):
*   egoAlterInt_i(j) = z1[s] * z2[j]         choice-only, candidate-dependent
*   tertius[j]       = mean_{k: A[k,j]>0} z[k]   both sub-models, candidate-ONLY
*                     = (z*A)[j] / indeg[j], with isolates (indeg[j]==0)
*                       imputed as the mean of tertius over non-isolates
*                       (goldfish's own documented imputation rule)
* both on the DEPENDENT network itself; default aggregateFun (mean) only.
*
* Real goldfish reference (Social_Evolution$calls, 84 actors, 439
* events, own `floor`/`gradeType` covariates - captured 2026-09-02,
* goldfish 1.6.12, see dev/dynam_unit10_tertius_egoalterint_crosscheck.R
* for the R side):
*   choice ~ tertius(floor):                     tertius=0.200163, logLik=-1892.8
*   choice ~ inertia + tertius(floor):            inertia=5.9087098, tertius=-0.0099198, logLik=-733.55
*   rate ~ tertius(floor):                        tertius=0.109599, logLik=-1932.6
*   choice ~ egoAlterInt(floor, gradeType):       egoAlterInt=0.025651, logLik=-1936.9
*   choice ~ inertia + alter(gradeType) + egoAlterInt(floor, gradeType):
*       inertia=5.683248, alter=0.829196, egoAlterInt=-0.176148, logLik=-711.11

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

preserve
import delimited "dev/dynam_crosscheck_actorcov.csv", clear varnames(1)
sort actor
assert _N == 84
mata: floorvec = st_data(1::84, "floor")'
* gradeType has real NAs in goldfish's own Social_Evolution$actors data,
* which R's write.csv renders as the literal string "NA" - forcing
* Stata's import delimited to read the whole column as a string. goldfish
* itself imputes missing attribute values with the MEAN of the
* non-missing values by default (estimate()'s own documented behavior:
* "For attributes missing values are replaced by the mean value") - so
* this crosscheck must match that exactly, not just destring-to-missing,
* to reproduce goldfish's own real reference numbers.
destring gradetype, replace force
summarize gradetype
replace gradetype = r(mean) if missing(gradetype)
mata: gradetypevec = st_data(1::84, "gradetype")'
restore

di as text "{hline 60}"
di as text "choice ~ tertius(floor)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0), zero8, zero8, zero8, zero8, floorvec, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b0", "V0", "ll0")
di as text "tertius=" %9.6f el(b0,1,1) "  logLik=" %9.4f ll0
scalar d0a = abs(el(b0,1,1) - 0.200163)
scalar d0b = abs(ll0 - (-1892.8))
di as text "abs diff: tertius=" %9.6f d0a "  logLik=" %9.4f d0b
assert d0a < 0.01
assert d0b < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + tertius(floor)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0), zero8, zero8, zero8, zero8, floorvec, zero8, zero8, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b1", "V1", "ll1")
di as text "inertia=" %9.6f el(b1,1,1) "  tertius=" %9.6f el(b1,1,2) "  logLik=" %9.4f ll1
scalar d1 = abs(el(b1,1,1) - 5.9087098)
scalar d2 = abs(el(b1,1,2) - (-0.0099198))
scalar d3 = abs(ll1 - (-733.55))
di as text "abs diff: inertia=" %9.6f d1 "  tertius=" %9.6f d2 "  logLik=" %9.4f d3
assert d1 < 0.01
assert d2 < 0.01
assert d3 < 0.1

di as text "{hline 60}"
di as text "rate ~ tertius(floor)"
di as text "{hline 60}"
mata: DynamRateFitMulti(eventmat, n_actors, (0,0,0,0,1), zero8, floorvec, 1e300, 1e300, 0, 0, 0, "b2", "V2", "ll2")
di as text "tertius=" %9.6f el(b2,1,1) "  logLik=" %9.4f ll2
scalar d4 = abs(el(b2,1,1) - 0.109599)
scalar d5 = abs(ll2 - (-1932.6))
di as text "abs diff: tertius=" %9.6f d4 "  logLik=" %9.4f d5
assert d4 < 0.01
assert d5 < 0.1

di as text "{hline 60}"
di as text "choice ~ egoAlterInt(floor, gradeType)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0), zero8, zero8, zero8, zero8, zero8, floorvec, gradetypevec, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b3", "V3", "ll3")
di as text "egoAlterInt=" %9.6f el(b3,1,1) "  logLik=" %9.4f ll3
scalar d6 = abs(el(b3,1,1) - 0.025651)
scalar d7 = abs(ll3 - (-1936.9))
di as text "abs diff: egoAlterInt=" %9.6f d6 "  logLik=" %9.4f d7
assert d6 < 0.01
assert d7 < 0.1

di as text "{hline 60}"
di as text "choice ~ inertia + alter(gradeType) + egoAlterInt(floor, gradeType)"
di as text "{hline 60}"
mata: DynamChoiceFitMulti(eventmat, n_actors, (1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0), zero8, zero8, zero8, gradetypevec, zero8, floorvec, gradetypevec, 1e300, 1e300, 0, 0, 0, 0, 0, J(0,0,0), J(0,0,0), "b4", "V4", "ll4")
assert colsof(b4) == 3
di as text "inertia=" %9.6f el(b4,1,1) "  alter=" %9.6f el(b4,1,2) "  egoAlterInt=" %9.6f el(b4,1,3) "  logLik=" %9.4f ll4
scalar d8 = abs(el(b4,1,1) - 5.683248)
scalar d9 = abs(el(b4,1,2) - 0.829196)
scalar d10 = abs(el(b4,1,3) - (-0.176148))
scalar d11 = abs(ll4 - (-711.11))
di as text "abs diff: inertia=" %9.6f d8 "  alter=" %9.6f d9 "  egoAlterInt=" %9.6f d10 "  logLik=" %9.4f d11
assert d8 < 0.01
assert d9 < 0.01
assert d10 < 0.01
assert d11 < 0.1

di as result "MATCH CONFIRMED (tertius under both sub-models, and egoAlterInt alone and alongside inertia+alter, all reproduce real goldfish)"
