cscript

* Time-varying covariates (option change()) against goldfish 1.6.12 (R 4.6): 42 random event logs with attribute change events for
* the choice, rate, rate-with-intercept and choice_coordination sub-models and every attribute effect (same, diff, sim, ego, alter,
* egoalterint, tertius), replace and increment changes, changes before the first and after the last event, changes at the time of
* an event (the event still sees the old value), several changes of one actor at one time, a window, several variables.
* The logs, the change lists and goldfish's estimates, standard errors and log likelihoods are the files of dev/tvc (made by
* dev/dynam_tvc_crosscheck.R). Tolerance: estimates, standard errors and ll to 1e-6 (relative; choice_coordination: estimates 1e-5,
* standard errors 10%, because its standard errors come from a numerical derivative). Then error paths and a no-change identity.
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
global TVDIR `"`c(pwd)'/dev/tvc"'
quietly do dev/dynam_tvc_crosscheck.do
runall

* ---- identity checks: a change event that sets an actor to its own value (or adds 0) gives the fit without change()
nwclear
import delimited using `"$TVDIR/nodes_c1.csv"', clear varnames(1) stringcols(1)
tempfile nodes
save `nodes'
import delimited using `"$TVDIR/events_c1.csv"', clear varnames(1) stringcols(2 3)
nwset sender receiver, eventtime(time) name(net) dropvars
rename _nwnode label
quietly merge 1:1 label using `nodes', keep(match master) nogen
rename label _nwnode
quietly nwdynam net, inertia alter(x)
matrix b0 = e(b)
scalar ll0 = e(ll)
* a change event that sets A1 to its own initial value changes nothing
local x1 = x[1]
preserve
clear
quietly set obs 1
quietly generate str3 node = "A1"
quietly generate str3 node2 = "A1"
quietly generate double time = 5.5
quietly generate double value = `x1'
quietly nwset node node2 value, eventtime(time) selfloop name(same_chg) dropvars
restore
quietly nwdynam net, inertia alter(x) change(x same_chg)
assert mreldif(e(b), b0) < 1e-12
assert reldif(e(ll), ll0) < 1e-12
di as result "PASS: a change that repeats the value leaves the fit unchanged"

* an increment of 0 also
preserve
clear
quietly set obs 1
quietly generate str3 node = "A1"
quietly generate str3 node2 = "A1"
quietly generate double time = 5.5
quietly generate double value = 0
quietly nwset node node2 value, eventtime(time) selfloop name(zero_chg) dropvars
restore
quietly nwdynam net, inertia alter(x) change(x zero_chg) increment
assert mreldif(e(b), b0) < 1e-12
di as result "PASS: an increment of 0 leaves the fit unchanged"

* option sequential runs with changes
quietly nwdynam net, inertia alter(x) change(x same_chg) sequential

* ---- error paths
capture nwdynam net, inertia alter(x) change(x)
assert _rc == 198
capture nwdynam net, inertia change(x same_chg)
assert _rc == 198
capture nwdynam net, inertia alter(x) change(y same_chg)
assert _rc == 111 | _rc == 198
* net itself holds events between different actors: not a change list
capture nwdynam net, inertia alter(x) change(x net)
assert _rc == 198
capture nwdynam net, inertia alter(x) increment
assert _rc == 198
capture nwdynam net, inertia alter(x) change(x same_chg x same_chg)
assert _rc == 198
* an actor that is not an actor of the fitted network
preserve
clear
quietly set obs 1
quietly generate str3 node = "Z9"
quietly generate str3 node2 = "Z9"
quietly generate double time = 5.5
quietly generate double value = 1
quietly nwset node node2 value, eventtime(time) selfloop name(bad_chg) dropvars
restore
capture nwdynam net, inertia alter(x) change(x bad_chg)
assert _rc == 198
di as result "PASS: error paths"
di as result "ALL OK"
