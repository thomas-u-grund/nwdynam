cscript

* Regression tests for two problems found while rewriting the book chapters.
* (1) "nwdynam version" as the FIRST command of a fresh session reported the native plugin as missing, because the
*     Mata library was not yet indexed (DynamNativePluginPath() "not found"). This test must therefore NOT load
*     unw_dynam.do: it runs against the built library lib/lnwdynam.mlib.
* (2) A fit whose effects are not identified (toy log of ten events) printed a table with standard errors in the
*     thousands and no warning; it now adds a note and returns e(unidentified).

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

* (1) first command of the session
nwdynam version
assert strpos(`"`r(plugin)'"', "missing") == 0
assert strpos(`"`r(plugin)'"', "dynam_sim") > 0
di as text "  PASS: nwdynam version finds the native plugin in a fresh session"

* (2) unidentified toy fit
nwclear
clear
input sender receiver t
1 2 1
1 3 2
2 1 3
1 2 4
3 2 5
2 3 6
1 3 7
3 1 8
2 1 9
1 3 10
end
nwset sender receiver, eventtime(t) directed name(chat)
nwdynam chat
assert strpos(" `e(unidentified)' ", " inertia ") > 0
assert strpos(" `e(unidentified)' ", " indeg ") > 0
assert strpos(" `e(unidentified)' ", " recip ") == 0
di as text "  PASS: unidentified effects are reported (`e(unidentified)')"

* an identified fit reports nothing
nwwebuse socialevolution, nwclear
nwset sender receiver, eventtime(t) directed name(calls)
nwdynam calls
assert `"`e(unidentified)'"' == ""
di as text "  PASS: identified fit has empty e(unidentified)"
di as result "ALL TESTS PASSED"
