cscript

* There is no Mata fallback: when the plugin file is missing nwdynam stops with r(601), when the installed plugin is too old for this
* version (or cannot be loaded) with r(3495); both with a message. With the plugin every sub-model fits.
* The test works in a scratch copy of the ado files (without lib/plugins) and puts a file of the plugin's name first on the adopath:
*  (a) no plugin anywhere              -> r(601), for every sub-model;
*  (b) a file that is not a plugin     -> r(3495);
*  (c) the plugin of this tree         -> fits for choice, rate, rate with intercept, choice_coordination.
* (The plugin of nwdynam 1.4.0 is tested in test_nwdynam_oldplugin.do: a plugin stays loaded for the rest of a Stata session.)
* Environment: NWCORE (core checkout, required).
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
local root `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
mata: mata mlib index

nwclear
clear
set obs 120
gen long sender = .
gen long receiver = .
gen eventtime = _n
mata:
mata set matastrict off
rseed(20261115)
n = 9
for (i=1; i<=120; i++) {
	pair = runiformint(1, 2, 1, n)
	while (pair[1] == pair[2]) pair = runiformint(1, 2, 1, n)
	st_store(i, "sender", pair[1])
	st_store(i, "receiver", pair[2])
}
end
tempfile evd
quietly save `evd'
nwset sender receiver, eventtime(eventtime) name(chk) dropvars
quietly use `evd', clear
nwset sender receiver, eventtime(eventtime) undirected name(chkc) dropvars
mata: st_local("pname", DynamNativePluginFilename())
mata: st_local("plugsrc", DynamNativePluginPath())
confirm file `"`plugsrc'"'

* scratch copy of the ado files, no plugin
tempfile __scx
local sc `"`__scx'_d"'
capture mkdir `"`sc'"'
capture erase `"`sc'/`pname'"'
foreach f in nwdynam.ado _nwdynam_libcheck.ado nwdynam_version.ado {
	quietly copy `"`root'/`f'"' `"`sc'/`f'"', replace
}
* remove the plugin directories from the adopath whatever slash style Stata stores them in (backslashes on Windows)
foreach __sub in macos unix windows {
	foreach __sep in "/" "\" {
		capture adopath - `"`root'`__sep'lib`__sep'plugins`__sep'`__sub'"'
	}
}
capture adopath - `"`root'"'
adopath ++ `"`sc'"'
capture program drop nwdynam
capture program drop _nwdynam_libcheck

* (a) no plugin
local i 0
foreach cmd in "nwdynam chk" "nwdynam chk, submodel(rate)" "nwdynam chk, submodel(rate) intercept" {
	local ++i
	capture noisily `cmd'
	assert _rc == 601
}
di as text "  PASS: no plugin file -> r(601) for choice, rate, rate with intercept"

* (b) a file that is no plugin
file open fh using `"`sc'/`pname'"', write text replace
file write fh "this is not a plugin" _n
file close fh
capture noisily nwdynam chk
assert _rc == 3495
capture noisily nwdynam chk, submodel(rate) intercept
assert _rc == 3495
di as text "  PASS: a file that cannot be loaded -> r(3495)"
erase `"`sc'/`pname'"'

* (c) the plugin of this tree
local sc2 `"`sc'/new"'
capture mkdir `"`sc2'"'
quietly copy `"`plugsrc'"' `"`sc2'/`pname'"', replace
capture shell chmod +x "`sc2'/`pname'"
adopath ++ `"`sc2'"'
nwdynam chk
nwdynam chk, submodel(rate)
nwdynam chk, submodel(rate) intercept
nwdynam chkc, submodel(choice_coordination)
di as text "  PASS: with the plugin all sub-models fit"
di "test_nwdynam_noplugin: ALL OK"
