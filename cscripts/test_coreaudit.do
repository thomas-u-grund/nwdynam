cscript

* Core 3.6.0 audit: no code change was needed in nwdynam; this test pins the behaviour on valued networks.
* - an event network made with a weight variable is valued; the weights are not used by the unweighted effects
*   (goldfish: the increment of an event is used only by the weighted* effects), so the estimates equal those of the
*   unweighted event network;
* - tie() reads whether a tie exists (goldfish: I(x_ij > 0)): a valued tie network gives the estimates of its
*   dichotomised network; negative values are not ties.

local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
quietly do unw_dynam.do

nwclear
clear
set obs 300
gen long sender = .
gen long receiver = .
gen eventtime = _n
gen double w = 1 + floor(runiform()*3)
mata:
mata set matastrict off
rseed(20260902)
for (i=1; i<=300; i++) {
	pair = runiformint(1, 2, 1, 10)
	while (pair[1] == pair[2]) pair = runiformint(1, 2, 1, 10)
	st_store(i, "sender", pair[1])
	st_store(i, "receiver", pair[2])
}
end
tempfile ev
save "`ev'"
nwset sender receiver, eventtime(eventtime) name(ev0)
use "`ev'", clear
nwset sender receiver w, eventtime(eventtime) name(ev1)
_nwsyntax ev1, max(1)
assert "`valued'" == "true"
nwdynam ev0, inertia recip
matrix B0 = e(b)
nwdynam ev1, inertia recip
assert mreldif(B0, e(b)) < 1e-3
di as text "  PASS: event weights are not used by the unweighted effects"

mata: rseed(7); M = (runiform(10,10):<.3); st_matrix("__S", M); st_matrix("__V", M :* (1:+floor(runiform(10,10):*4)))
nwset, mat(__S) directed name(st0) nooutput
nwset, mat(__V) directed name(stv) nooutput
_nwsyntax stv, max(1)
assert "`valued'" == "true"
nwdynam ev0, inertia tie(st0)
matrix B1 = e(b)
nwdynam ev0, inertia tie(stv)
assert mreldif(B1, e(b)) < 1e-5
* negative values are not ties: -V gives the same as a network without those ties
mata: Mn = st_matrix("__V"); Mn = Mn :* (2*(runiform(10,10):<.5) :- 1); st_matrix("__Neg", Mn); st_matrix("__Pos", (Mn :> 0))
nwset, mat(__Neg) directed name(stn) nooutput
nwset, mat(__Pos) directed name(stp) nooutput
nwdynam ev0, inertia tie(stn)
matrix B2 = e(b)
nwdynam ev0, inertia tie(stp)
assert mreldif(B2, e(b)) < 1e-5
di as text "  PASS: tie() of a valued network equals tie() of its dichotomised network"
di as text "test_coreaudit (nwdynam) PASS"
