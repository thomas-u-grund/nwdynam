cscript
local nwcore : env NWCORE
adopath ++ `"`nwcore'"'
capture adopath ++ `"`nwcore'/lib"'
adopath ++ `"`c(pwd)'"'
adopath ++ `"`c(pwd)'/lib"'
adopath ++ `"`c(pwd)'/lib/plugins/macos"'
global COMPDIR `"`c(pwd)'/dev/comp"'
quietly do dev/dynam_comp_crosscheck.do
runall
