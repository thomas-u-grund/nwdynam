* Platform detection: console Stata on a Mac (for example StataNow 19.5) reports c(os) = "Unix" but
* c(machine_type) starts with "Mac"; the plugin lookup must then ask for the macos file, never the unix (Linux) one.
* Run from the root of the repository.
clear all
adopath ++ "`c(pwd)'"
adopath ++ "`c(pwd)'/lib"
mata: mata mlib index
mata: st_local("fn", DynamNativePluginFilename()); st_local("sd", DynamNativePluginSubdir())
local want "unix"
if c(os) == "Windows" local want "windows"
if c(os) == "MacOSX" | substr(c(machine_type), 1, 3) == "Mac" local want "macos"
di "c(os) = " c(os) ", c(machine_type) = " c(machine_type) ", plugin file = `fn', subdirectory = `sd'"
assert "`sd'" == "`want'"
assert "`fn'" == "nwdynam_sim_`want'.plugin" | "`fn'" == "dynam_sim_`want'.plugin"
di "PASS: platform detection"
