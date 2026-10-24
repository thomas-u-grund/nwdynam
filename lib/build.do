* Builds lib/lnwdynam.mlib from the Mata sources of this package. Run from the root of the repository:
*   StataBE -e do lib/build.do
* The library contains only the functions of this package (the core is not loaded while compiling).
version 12
clear all
capture mata: mata clear
do unw_dynam.do
local root `"`c(pwd)'"'
mata: mata mlib create lnwdynam, dir(`"`root'/lib"') replace
mata: mata mlib add lnwdynam *(), dir(`"`root'/lib"')
mata: mata mlib index
