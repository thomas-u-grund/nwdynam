global WINDIR "`c(pwd)'/dev/win"
quietly do dev/dynam_win_crosscheck.do
quietly runcase f11 choice "inertia tertius(y) four alter(x) indeg" two "x y" 0 0
matrix list e(b)
quietly runcase f12 choice "inertia same(g) diff(x) sim(x) alter(x) egoalterint(x y) tertius(x) trans" one "" 0 0
matrix list e(b)
