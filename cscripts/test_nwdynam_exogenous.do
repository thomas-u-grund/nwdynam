cscript

* Effects on other networks of nwdynam against goldfish 1.6.12 (R 4.6): (a) tie() of static and changing networks, weighted, windowed, ignoreRep,
* (b) trans/cycle/commonSender/commonReceiver through exogenous networks, (c) mixed closures (mixedTrans, mixedCycle, mixedCommonSender,
* mixedCommonReceiver), (d) tertiusDiff, (e) ignoreRep on the effects of the dependent network, (f) choice_coordination.
* The logs and goldfish's estimates, standard errors and log likelihoods are the files of dev/ex (made by dev/dynam_ex_crosscheck.R, regenerate with
* `Rscript dev/dynam_ex_crosscheck.R dev/ex`). Tolerance: estimates, standard errors and ll to 1e-6 (relative); choice_coordination: estimates 5e-4,
* standard errors 10% (flat likelihood, numerical derivative). Then fixed reference values and the error paths of what stays refused.
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
global WINDIR `"`c(pwd)'/dev/ex"'
global WINSEC ""
global WINONLY ""
local sec : env EXSEC
if `"`sec'"' != "" global WINSEC `"`sec'"'
local only : env EXONLY
if `"`only'"' != "" global WINONLY `"`only'"'
quietly do dev/dynam_ex_crosscheck.do
runall
