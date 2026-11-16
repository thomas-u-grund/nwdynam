global WINDIR "`c(pwd)'/dev/win"
quietly do dev/dynam_win_crosscheck.do
foreach o in "inertiawindow(2) weightedinertia" "inertiawindow(2)" "indegwindow(2) inertiawindow(2) weightedinertia" {
	quietly runcase c37 choice_coordination "`o'" und "" 0 0
	di "RES `o' : " e(b)[1,1] " " e(b)[1,2] " ll " e(ll)
}
