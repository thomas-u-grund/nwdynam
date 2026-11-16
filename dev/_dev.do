local sec : env WSEC
local only : env WONLY
global WINDIR "`c(pwd)'/dev/win"
global WINSEC "`sec'"
global WINONLY "`only'"
quietly do dev/dynam_win_crosscheck.do
runall
