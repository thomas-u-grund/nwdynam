* Builds lib/lnwdynam.mlib from the Mata sources of this package. Run from the root of the repository:
*   StataBE -e do lib/build.do
* The library contains only the functions of this package (the core is not loaded while compiling).
version 12
* BUILD GUARD: Mata libraries must be built with the OLDEST supported Stata (currently Stata 17), never with a newer one: a library
* compiled by a newer Stata is silently ignored by older Stata versions ("compiled by Stata 18.0, is too new to be run by this version
* of Stata and so was ignored"), and every function is then "not found" (r(3499)). Set NW_ALLOW_NEWER_STATA=1 only for experiments
* whose library is never released.
local nw_minbuild 17
local nw_allow : env NW_ALLOW_NEWER_STATA
if c(stata_version) > `nw_minbuild' & "`nw_allow'" != "1" {
    di as error "this library must be built with Stata `nw_minbuild' (the oldest supported version), but this is Stata " c(stata_version) ". A library built by a newer Stata cannot be used by older Stata versions."
    di as error "Run the build with Stata `nw_minbuild' (see RELEASING.md, 'Build the Mata libraries with the oldest supported Stata'); NW_ALLOW_NEWER_STATA=1 overrides this for experiments."
    exit 198
}
clear all
capture mata: mata clear
do unw_dynam.do
local root `"`c(pwd)'"'
mata: mata mlib create lnwdynam, dir(`"`root'/lib"') replace
mata: mata mlib add lnwdynam *(), dir(`"`root'/lib"')
mata: mata mlib index
