cscript

* The plugins of nwdynam 1.4.0 (protocol 1) and 1.6.0 (protocol 3, before composition change), taken from the git history (skipped when
* git or that commit is not available), are too old for this version: every sub-model stops with r(3495) and a message.
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
nwset sender receiver, eventtime(eventtime) name(chk)
quietly use `evd', clear
nwset sender receiver, eventtime(eventtime) undirected name(chkc)
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
capture adopath - `"`root'/lib/plugins/macos"'
capture adopath - `"`root'/lib/plugins/unix"'
capture adopath - `"`root'/lib/plugins/windows"'
capture adopath - `"`root'"'
adopath ++ `"`sc'"'
capture program drop nwdynam
capture program drop _nwdynam_libcheck

if "`c(os)'" == "MacOSX" {
	foreach __rev in ae51be1 4f3188a {
		quietly shell git -C "`root'" show `__rev':lib/plugins/macos/`pname' > "`sc'/`pname'" 2>/dev/null
		capture confirm file `"`sc'/`pname'"'
		local __sz 0
		if _rc == 0 {
			mata: __fh = fopen(st_local("sc") + "/" + st_local("pname"), "r"); fseek(__fh, 0, 1); st_local("__sz", strofreal(ftell(__fh))); fclose(__fh)
			capture mata: mata drop __fh
		}
		if `__sz' > 0 {
			adopath ++ `"`sc'"'
			quietly shell chmod +x "`sc'/`pname'"
			foreach cmd in "nwdynam chk" "nwdynam chk, submodel(rate)" "nwdynam chk, submodel(rate) intercept" {
				capture noisily `cmd'
				assert _rc == 3495
			}
			di as text "  PASS: the plugin of commit `__rev' -> r(3495) for choice, rate, rate with intercept"
			erase `"`sc'/`pname'"'
		}
		else di as text "  (git history not available; the check of commit `__rev' is skipped)"
	}
}

di "test_nwdynam_oldplugin: ALL OK"
