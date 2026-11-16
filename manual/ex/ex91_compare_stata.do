* Chapter 11: the nwdynam side of the comparison with goldfish: fits the same models and writes ex/r/nwdynam_results.csv
do setup.doh
tempname post
tempfile res
postfile `post' str20 model str20 term double(estimate se loglik) using `res', replace
capture program drop store
program store
    args post tag
    matrix b = e(b)
    matrix V = e(V)
    local names : colnames b
    local j = 0
    foreach nm of local names {
        local ++j
        post `post' ("`tag'") ("`nm'") (b[1,`j']) (sqrt(V[`j',`j'])) (e(ll))
    }
end
*# cmp_soc
import delimited using "r/socevo.csv", clear
quietly levelsof sender, local(a)
quietly levelsof receiver, local(b)
local seen : list a | b
display "actors in the log: " `: word count `seen''
tostring sender receiver, replace
nwset sender receiver, directed eventtime(t) name(soc)
local miss
forvalues i = 1/84 {
    if !`: list posof "`i'" in seen' local miss `miss',`i'
}
if "`miss'" != "" {
    local miss = substr("`miss'", 2, .)
    nwaddnodes soc, nodenames(`miss')
}
nwset
nwdynam soc, inertia recip indeg
store `post' soc_choice
nwdynam soc, submodel(rate) indeg outdeg
store `post' soc_rate
nwdynam soc, submodel(rate) intercept indeg outdeg
store `post' soc_rate_int
*# cmp_email
mkmail 5
nwdynam mail, inertia recip indeg
store `post' em_choice
nwdynam mail, inertia recip indeg trans
store `post' em_choice_trans
nwdynam mail, inertiawindow(0.5) recip indeg
store `post' em_choice_win
nwdynam mail, submodel(rate) indeg outdeg
store `post' em_rate
nwdynam mail, submodel(rate) intercept indeg outdeg
store `post' em_rate_int
postclose `post'
use `res', clear
export delimited using "r/nwdynam_results.csv", replace
*# end
