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

* ---- fixed reference values (goldfish 1.6.12)
quietly runcase a01 choice "netterms(tie(x1))" one "x1=dyndir"
assert reldif(el(e(b), 1, 1), -0.410100241262744) < 1e-6 & reldif(e(ll), -89.0111350828585) < 1e-6
quietly runcase a02 choice "netterms(tie(x3, weighted))" one "x3=dyndir"
assert reldif(el(e(b), 1, 1), 0.246327826406341) < 1e-6 & reldif(e(ll), -86.4071361404911) < 1e-6
quietly runcase a03 choice "netterms(tie(x1, window(6)))" one "x1=dyndir"
assert reldif(el(e(b), 1, 1), -0.344120317923601) < 1e-6 & reldif(e(ll), -89.4777219075591) < 1e-6
quietly runcase a04 choice "netterms(tie(x2))" one "x2=statdir"
assert reldif(el(e(b), 1, 1), -0.0307326589319568) < 1e-6 & reldif(e(ll), -89.5828581105147) < 1e-6
quietly runcase b06 choice "netterms(four(x1))" one "x1=dyndir"
assert reldif(el(e(b), 1, 1), -0.301493583746457) < 1e-6 & reldif(e(ll), -87.739426102554) < 1e-6
quietly runcase b07 choice "netterms(four(x2))" one "x2=statdir"
assert reldif(el(e(b), 1, 1), -0.183546229954431) < 1e-6 & reldif(e(ll), -89.4669720626772) < 1e-6
quietly runcase c01 choice "netterms(mixedtrans(net x1))" one "x1=dyndir"
assert reldif(el(e(b), 1, 1), -0.194811925434428) < 1e-6 & reldif(e(ll), -89.2492344853782) < 1e-6
quietly runcase d03 choice "netterms(tertiusdiff(y x2))" one "x2=statdir"
assert reldif(el(e(b), 1, 1), -0.0724097522508217) < 1e-6 & reldif(e(ll), -89.5596725611363) < 1e-6
quietly runcase d05 choice "inertia netterms(tertiusdiff(x))" one ""
assert reldif(el(e(b), 1, 1), 0.37293016030384) < 1e-6 & reldif(el(e(b), 1, 2), -0.349414526558267) < 1e-6 & reldif(e(ll), -114.636777189654) < 1e-6
quietly runcase e01 choice "inertia indeg ignorerep(indeg)" one ""
assert reldif(el(e(b), 1, 1), 0.389400900792924) < 1e-6 & reldif(el(e(b), 1, 2), -0.127217159443728) < 1e-6 & reldif(e(ll), -87.6319810027618) < 1e-6
quietly runcase e02 choice "inertia trans ignorerep(trans)" one ""
assert reldif(el(e(b), 1, 1), 0.795832743777731) < 1e-6 & reldif(el(e(b), 1, 2), 0.152468007824017) < 1e-6 & reldif(e(ll), -87.79730071891) < 1e-6
quietly runcase f04 choice_coordination "netterms(mixedtrans(net x1))" und "x1=dynund"
assert reldif(el(e(b), 1, 1), -0.542253774979325) < 5e-4 & reldif(e(ll), -145.536414440647) < 5e-4
quietly runcase f08 choice_coordination "inertia four ignorerep(four)" und ""
assert reldif(el(e(b), 1, 1), 0.483669174778511) < 5e-4 & reldif(el(e(b), 1, 2), 0.0219127766169733) < 5e-4 & reldif(e(ll), -197.235977233724) < 5e-4
di as result "PASS: fixed reference values"

* ---- tertiusDiff: where goldfish stops with an error (a candidate loses its last in-neighbour: window 1, 2, 3 on this log), nwdynam stops too;
* with a window longer than the log nothing expires and it agrees (goldfish 1.6.12: -0.583924986997476, ll -88.7209609019643)
foreach w in 1 2 3 {
	capture runcase d01 choice "netterms(tertiusdiff(x, window(`w')))" one ""
	if _rc == 0 {
		di as error "tertiusdiff window(`w'): expected an error"
		exit 9
	}
}
quietly runcase d01 choice "netterms(tertiusdiff(x, window(100)))" one ""
assert reldif(el(e(b), 1, 1), -0.583924986997476) < 1e-6 & reldif(e(ll), -88.7209609019643) < 1e-6
di as result "PASS: tertiusdiff windows"

* ---- error paths: what goldfish refuses or does not give stays refused (all r(198))
foreach spec in ///
	"a01|choice|netterms(foo(x1))|x1=dyndir" ///
	"a01|choice|netterms(tie(x1, bogus))|x1=dyndir" ///
	"a01|choice|netterms(tie(x1 x1))|x1=dyndir" ///
	"a01|choice|netterms(tie)|x1=dyndir" ///
	"a01|choice|netterms(mixedtrans(x1))|x1=dyndir" ///
	"a01|choice|netterms(mixedtrans(x1 x1))|x1=dyndir" ///
	"a01|choice|netterms(trans(x1, weighted))|x1=dyndir" ///
	"a01|choice|netterms(tie(x1, window(0)))|x1=dyndir" ///
	"a01|choice|netterms(tie(x1, window(-3)))|x1=dyndir" ///
	"a01|choice|netterms(tertiusdiff(x x1 x1))|x1=dyndir" ///
	"a01|choice|netterms(tie(x1) tie(x1))|x1=dyndir" ///
	"a01|rate|netterms(tie(x1))|x1=dyndir" ///
	"a01|choice|ignorerep(same)|x1=dyndir" ///
	"a01|choice|inertia ignorerep(trans)|x1=dyndir" ///
	"a01|choice|ignorerep(nosuchthing)|x1=dyndir" ///
	"a01|rate|ignorerep(indeg)|x1=dyndir" ///
	"a01|choice|tertius(x) netterms(tertiusdiff(x)) change(x nochg)|x1=dyndir" ///
	"f01|choice_coordination|netterms(cycle(x1))|x1=dynund" ///
	"f01|choice_coordination|netterms(mixedcycle(net x1))|x1=dynund" ///
	"f01|choice_coordination|netterms(commonsender(x1))|x1=dynund" {
	tokenize `"`spec'"', parse("|")
	local kd one
	if "`3'" == "choice_coordination" local kd und
	capture runcase `1' `3' "`5'" `kd' "`7'"
	local rc = _rc
	if `rc' != 198 {
		di as error "`spec': expected r(198), got r(`rc')"
		exit 9
	}
}
* a network that is not there, a static network with another number of actors, an event network with an actor the fitted network does not have
quietly runcase a01 choice "netterms(tie(x1))" one "x1=dyndir"
capture nwdynam net, netterms(tie(nosuchnetwork))
assert _rc != 0
nwrandom 5, prob(.3) name(smallnet) directed
capture nwdynam net, netterms(tie(smallnet))
assert _rc == 198
clear
input str3 s str3 r double t
"A1" "A2" 1
"A2" "A3" 2
"A3" "A1" 3
"A1" "A3" 4
"A2" "A1" 5
end
nwset s r, eventtime(t) directed name(tinynet)
clear
input str3 s str3 r double t
"A1" "Z9" 1
"A2" "A3" 2
end
nwset s r, eventtime(t) directed name(strangenet)
capture nwdynam tinynet, netterms(tie(strangenet))
assert _rc == 198
* two-mode networks
clear
input str3 s str3 r double t
"P1" "O1" 1
"P2" "O1" 2
"P1" "O2" 3
"P2" "O2" 4
"P3" "O1" 5
"P3" "O2" 6
end
nwset s r, twomode eventtime(t) name(twonet)
capture nwdynam twonet, netterms(tie(tinynet))
assert _rc == 198
di as result "PASS: error paths (what stays refused)"
di as result "ALL OK"
