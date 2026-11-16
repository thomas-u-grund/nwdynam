cscript

* Composition change (options present() and compchange()) against goldfish 1.6.12 (R 4.6): actors who join or leave the network.
* 1. 67 random and fixed event logs for the choice, rate, rate-with-intercept and choice_coordination sub-models and one- and two-mode
*    networks: joins and leaves before the first and after the last event, at the time of an event, several changes of one actor at
*    one time, actors absent during the whole log, redundant changes, windows, time-varying covariates together with composition
*    change. The logs, the change lists and goldfish's estimates, standard errors and log likelihoods are the files of dev/comp (made by
*    dev/dynam_comp_crosscheck.R). Tolerance: estimates, standard errors and ll to 1e-6 (relative; choice_coordination: estimates
*    1e-5, standard errors 10%, because its standard errors come from a numerical derivative).
* 2. a fixed reference (goldfish values typed in) with opportunities() and a join/leave network, 3. identities, 4. error paths.
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
capture adopath ++ `"`nwcore'/lib"'
adopath ++ `"`c(pwd)'"'
adopath ++ `"`c(pwd)'/lib"'
adopath ++ `"`c(pwd)'/lib/plugins/macos"'
global COMPDIR `"`c(pwd)'/dev/comp"'
quietly do dev/dynam_comp_crosscheck.do
runall

* ---- 2. fixed reference: choice ~ inertia + indeg, 8 actors, 40 events; A7 is absent until its join, A8 leaves and returns
*         (the absence of A7 at the start is a leave event before the first event); goldfish 1.6.12 values (estimates and ll):
*         without opportunities 0.09813994 0.1138617 ll -73.33197; with a per-event opportunity list -0.04669544 0.1190584 ll -67.61801
nwclear
import delimited using `"$COMPDIR/comp_opp.csv"', clear varnames(1) stringcols(2)
quietly set obs `=_N + 1'
quietly replace time = 0.01 in l
quietly replace node = "A7" in l
quietly replace value = 0 in l
generate node2 = node
quietly nwset node node2 value, eventtime(time) selfloop name(comp)
import delimited using `"$COMPDIR/events_opp.csv"', clear varnames(1) stringcols(2 3)
quietly nwset sender receiver, eventtime(time) directed name(net)
foreach a in A1 A2 A3 A4 A5 A6 A7 A8 {
	capture quietly nwaddnodes net, nodenames(`a')
}
quietly nwdynam net, inertia indeg compchange(comp)
assert abs(el(e(b), 1, 1) - 0.09813994) < 1e-5 & abs(el(e(b), 1, 2) - 0.1138617) < 1e-5 & abs(e(ll) - (-73.33197)) < 1e-4
assert e(ncompchange) == 4
assert "`e(compchange)'" == "comp"
di as result "PASS: fixed reference without opportunities"
import delimited using `"$COMPDIR/avail_opp.csv"', clear varnames(1)
quietly nwdynam net, inertia indeg compchange(comp) opportunities(event actor)
assert abs(el(e(b), 1, 1) - (-0.04669544)) < 1e-5 & abs(el(e(b), 1, 2) - 0.1190584) < 1e-5 & abs(e(ll) - (-67.61801)) < 1e-4
di as result "PASS: fixed reference with opportunities() (the risk set is the intersection of the opportunities and the present actors)"

* ---- 3. identities: everybody present (present() all 1, or joins of actors who are present) equals the fit without composition change
nwclear
import delimited using `"$COMPDIR/nodes_c01.csv"', clear varnames(1) stringcols(1)
tempfile ndf
save `ndf'
import delimited using `"$COMPDIR/events_c01.csv"', clear varnames(1) stringcols(2 3)
quietly nwset sender receiver, eventtime(time) directed name(net)
foreach a in A1 A2 A3 A4 A5 A6 {
	capture quietly nwaddnodes net, nodenames(`a')
}
quietly nwdynam net, inertia indeg
matrix b0 = e(b)
scalar ll0 = e(ll)
quietly nwdynam net, inertia indeg submodel(choice)
assert mreldif(e(b), b0) == 0
rename _nwnode label
quietly merge 1:1 label using `ndf', nogen
rename label _nwnode
generate byte allpres = 1
quietly nwdynam net, inertia indeg present(allpres)
assert mreldif(e(b), b0) < 1e-12 & reldif(e(ll), ll0) < 1e-12
di as result "PASS: present() = 1 for everybody leaves the choice fit unchanged"
* joins of present actors, and a leave and a return that no event of the actor spans... (redundant changes)
preserve
clear
quietly set obs 3
generate str3 node = "A1"
generate str3 node2 = "A1"
generate double time = 3
generate double value = 1
quietly replace time = 1e6 in 2
quietly replace value = 0 in 2
quietly replace time = 2e6 in 3
quietly replace value = 1 in 3
quietly nwset node node2 value, eventtime(time) selfloop name(redund)
restore
quietly nwdynam net, inertia indeg compchange(redund)
assert mreldif(e(b), b0) < 1e-12 & reldif(e(ll), ll0) < 1e-12
di as result "PASS: redundant joins and a leave after the last event leave the fit unchanged"
quietly nwdynam net, submodel(rate) indeg outdeg
matrix b1 = e(b)
scalar ll1 = e(ll)
quietly nwdynam net, submodel(rate) indeg outdeg compchange(redund) present(allpres)
assert mreldif(e(b), b1) < 1e-12 & reldif(e(ll), ll1) < 1e-12
quietly nwdynam net, submodel(rate) indeg outdeg intercept
matrix b2 = e(b)
scalar ll2 = e(ll)
quietly nwdynam net, submodel(rate) indeg outdeg intercept compchange(redund) present(allpres)
assert mreldif(e(b), b2) < 1e-12 & reldif(e(ll), ll2) < 1e-12
di as result "PASS: the same for the rate sub-model with and without intercept"
* A1 has events: absent from the start is an error (section 4)
generate byte abs1 = 1
quietly replace abs1 = 0 if _nwnode == "A1"

* ---- 4. error paths
* an event of an absent actor (goldfish: "the node is not present at time"), at the start, after a leave, and with the leave at the time of the event
capture noisily nwdynam net, inertia indeg present(abs1)
assert _rc == 198
di as result "PASS: an actor absent at the start who has an event -> r(198)"
preserve
clear
quietly set obs 1
generate str3 node = "A1"
generate str3 node2 = "A1"
generate double time = 0.5
generate double value = 0
quietly nwset node node2 value, eventtime(time) selfloop name(leave1)
restore
capture noisily nwdynam net, inertia indeg compchange(leave1)
assert _rc == 198
di as result "PASS: a leave before the actor's events -> r(198)"
* (times that must equal the event times exactly are float, as the times of the imported events)
* the join at the time of the first event of the actor makes it present for that event (goldfish: changes at the time of an event count)
preserve
import delimited using `"$COMPDIR/events_c01.csv"', clear varnames(1) stringcols(2 3)
quietly keep if sender == "A1" | receiver == "A1"
local tfirst = time[1]
local tlast = time[_N]
clear
quietly set obs 2
generate str3 node = "A1"
generate str3 node2 = "A1"
generate float time = 0.5
generate double value = 0
quietly replace time = `tfirst' in 2
quietly replace value = 1 in 2
quietly nwset node node2 value, eventtime(time) selfloop name(joinatfirst)
restore
quietly nwdynam net, inertia indeg compchange(joinatfirst)
di as result "PASS: a join at the time of the actor's first event is allowed"
* ... but a leave at the time of the actor's last event is not
preserve
clear
quietly set obs 1
generate str3 node = "A1"
generate str3 node2 = "A1"
generate float time = `tlast'
generate double value = 0
quietly nwset node node2 value, eventtime(time) selfloop name(leaveatlast)
restore
capture noisily nwdynam net, inertia indeg compchange(leaveatlast)
assert _rc == 198
di as result "PASS: a leave at the time of the actor's event -> r(198)"
* malformed change networks
capture program drop mkchg
program mkchg
	args nm n1 n2 t v
	preserve
	clear
	quietly set obs 1
	quietly generate str3 node = "`n1'"
	quietly generate str3 node2 = "`n2'"
	quietly generate double time = `t'
	quietly generate double value = `v'
	if "`n1'" == "`n2'" quietly nwset node node2 value, eventtime(time) selfloop name(`nm')
	else quietly nwset node node2 value, eventtime(time) name(`nm')
	restore
end
mkchg notself A1 A2 3 1
mkchg self1 A1 A1 3 1
mkchg val2 A1 A1 3 2
mkchg unknown Z9 Z9 3 1
foreach bad in notself val2 unknown {
	capture noisily nwdynam net, inertia indeg compchange(`bad')
	assert _rc == 198
}
di as result "PASS: sender != receiver, a value other than 0/1, an unknown actor in compchange() -> r(198)"
capture noisily nwdynam net, inertia indeg compchange(self1 val2)
assert _rc == 198
preserve
clear
quietly set obs 2
generate str3 node = "A1"
generate str3 node2 = "A2"
quietly nwset node node2, name(staticnet)
restore
capture noisily nwdynam net, inertia indeg compchange(staticnet)
assert _rc == 198
di as result "PASS: compchange() of two names or of a network without event times -> r(198)"
* present(): values other than 0/1 or missing
generate double badpres = 2
capture noisily nwdynam net, inertia indeg present(badpres)
assert _rc == 198
quietly replace badpres = .
capture noisily nwdynam net, inertia indeg present(badpres)
assert _rc == 198
di as result "PASS: present() with a value other than 0/1, or a missing value -> r(198)"
* a covariate change event (change()) of an absent actor is an error, as in goldfish; at the time of the join it is allowed
preserve
import delimited using `"$COMPDIR/events_c01.csv"', clear varnames(1) stringcols(2 3)
quietly keep if sender == "A2" | receiver == "A2"
local t2 = time[1]
clear
quietly set obs 2
generate str3 node = "A2"
generate str3 node2 = "A2"
generate float time = 0.1
generate double value = 0
quietly replace time = `t2' in 2
quietly replace value = 1 in 2
quietly nwset node node2 value, eventtime(time) selfloop name(a2comp)
clear
quietly set obs 1
generate str3 node = "A2"
generate str3 node2 = "A2"
generate float time = `t2' / 2
generate double value = 1.5
quietly nwset node node2 value, eventtime(time) selfloop name(a2chgbad)
clear
quietly set obs 1
generate str3 node = "A2"
generate str3 node2 = "A2"
generate float time = `t2'
generate double value = 1.5
quietly nwset node node2 value, eventtime(time) selfloop name(a2chgok)
restore
capture noisily nwdynam net, inertia alter(x) change(x a2chgbad) compchange(a2comp)
assert _rc == 198
quietly nwdynam net, inertia alter(x) change(x a2chgok) compchange(a2comp)
di as result "PASS: a covariate change of an absent actor -> r(198); at the time of its join it is allowed"
di "test_nwdynam_composition: ALL OK"
