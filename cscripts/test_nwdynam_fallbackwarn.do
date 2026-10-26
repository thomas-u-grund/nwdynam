cscript

* nwdynam prints one warning line when the compiled plugin is missing and the original structural effect set (choice
* and rate sub-models) runs on the Mata fallback; no warning with the plugin; no warning for Mata-only effect sets.
* Run this do-file twice: once in a tree with lib/plugins (expects no warning, writes the reference coefficients) and once in
* a copy of the tree without the plugin files (expects exactly one warning per call, coefficients equal to the reference).
* Environment: NWCORE (core checkout, required), NWDYN_REF (path of the reference file, optional).

local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
adopath ++ `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
local ref : env NWDYN_REF

nwclear
clear
set obs 300
gen long sender = .
gen long receiver = .
gen eventtime = _n
mata:
mata set matastrict off
rseed(20261026)
n = 12
for (i=1; i<=300; i++) {
	pair = runiformint(1, 2, 1, n)
	while (pair[1] == pair[2]) pair = runiformint(1, 2, 1, n)
	st_store(i, "sender", pair[1])
	st_store(i, "receiver", pair[2])
}
end
nwset sender receiver, eventtime(eventtime) name(chatlog)

mata: st_local("avail", strofreal(DynamNativeAvailable()))
di as text "plugin available: `avail'"

tempfile lg
program define _countwarn, rclass
	args file
	tempname fh
	local n = 0
	file open `fh' using `"`file'"', read text
	file read `fh' line
	while r(eof) == 0 {
		if strpos(`"`macval(line)'"', "Warning: the compiled plugin is not used") local ++n
		file read `fh' line
	}
	file close `fh'
	return scalar n = `n'
end

* choice (default), rate: eligible engines
foreach sm in choice rate {
	capture log close wl
	log using `"`lg'_`sm'.log"', text replace name(wl)
	nwdynam chatlog, submodel(`sm')
	log close wl
	_countwarn `"`lg'_`sm'.log"'
	local nw_`sm' = r(n)
	tempname b_`sm'
	matrix `b_`sm'' = e(b)
	di as text "  `sm': warnings = `nw_`sm''"
	assert `nw_`sm'' == cond(`avail', 0, 1)
}

* Mata-only effect set (inertia only is a genuine subset): never a warning
capture log close wl
log using `"`lg'_sub.log"', text replace name(wl)
nwdynam chatlog, inertia
log close wl
_countwarn `"`lg'_sub.log"'
assert r(n) == 0
di as text "  subset engine: no warning (Mata-only by design)"

* coefficients identical to the plug-in run
mata: bc = st_matrix(st_local("b_choice")); br = st_matrix(st_local("b_rate"))
if `"`ref'"' != "" {
	capture confirm file `"`ref'"'
	if _rc {
		mata: fh = fopen(st_local("ref"), "w"); fput(fh, invtokens(strofreal(bc, "%21.17g"))); fput(fh, invtokens(strofreal(br, "%21.17g"))); fclose(fh)
		di as text "  reference coefficients written"
	}
	else {
		mata: fh = fopen(st_local("ref"), "r"); rc = strtoreal(tokens(fget(fh))); rr = strtoreal(tokens(fget(fh))); fclose(fh)
		mata: assert(max(abs(bc :- rc)) < 1e-6); assert(max(abs(br :- rr)) < 1e-6)
		di as text "  coefficients equal to the plug-in run (1e-6)"
	}
}
di as result "ALL OK test_nwdynam_fallbackwarn"
