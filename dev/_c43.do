global WINDIR "`c(pwd)'/dev/win"
quietly do dev/dynam_win_crosscheck.do
import delimited using "$WINDIR/results_c.csv", clear varnames(1)
list if case=="c43"
list if case=="c27"
foreach c in c43 c27 {
	import delimited using `"$WINDIR/cases_c.csv"', clear varnames(1) stringcols(1 2 3 4 5)
	levelsof opts if id=="`c'", local(o) clean
	di "`o'"
	forvalues s = 1/3 {
		set seed `s'
		runcase `c' choice_coordination "`o' seed(`s')" und "" 0 0
		matrix list e(b)
		di e(ll)
	}
}
