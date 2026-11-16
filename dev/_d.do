global WINDIR "`c(pwd)'/dev/win"
quietly do dev/dynam_win_crosscheck.do
runcase b92 rate "tertius(x) tertiuswindow(2) intercept" one "x" 0 1
