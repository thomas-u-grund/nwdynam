* Chapters 7 and 8: effects on other networks, time-varying covariates, composition change
do setup.doh
* who is in the log before day 4.5 (used for the composition example below)
quietly nwwebuse manufacturing, nwclear
keep if t < 4.5
keep sender receiver
gen label = sender
tempfile s1
quietly save `s1'
replace label = receiver
append using `s1'
keep label
quietly duplicates drop
gen byte seen = 1
tempfile early
quietly save `early'
quietly nwwebuse manufacturing_nodes, nwclear
encode division, generate(div)
tempfile actors
quietly save `actors'
mkmail 5
*# sp_attr
nwload mail, xvars
gen long row = _n
gen label = _nwnode
quietly merge m:1 label using `actors', keep(master match) nogenerate
sort row
list label supervisor level div in 1/4
*# sp_hier
* the reporting hierarchy as an ordinary (static) network in the actor order of mail: the tie i -> j is present when j is the supervisor of i
mata:
lab = st_sdata(., "label"); sup = st_sdata(., "supervisor"); n = rows(lab); M = J(n, n, 0)
for (i=1; i<=n; i++) {
    j = selectindex(lab :== sup[i])
    if (length(j)) M[i, j[1]] = 1
}
st_matrix("HM", M)
end
nwset, mat(HM) name(hier) directed
nwset
*# sp_tie
nwdynam mail, inertia recip indeg netterms(tie(hier))
*# sp_mixed
nwdynam mail, inertia recip indeg netterms(tie(hier) mixedtrans(hier mail))
*# sp_tdiff
nwdynam mail, inertia recip indeg netterms(tertiusdiff(level))
*# sp_ignorerep
nwdynam mail, indeg trans ignorerep(trans)
*# sp_changedata
* artificial change events, for illustration: ten employees, chosen at random, move one level up on day 4.2
set seed 2026
gen double u = runiform()
sort u
gen byte moved = _n <= 10
sort row
preserve
keep if moved
gen newlevel = level - 1
gen label2 = label
gen double tchg = 4.2
nwset label label2 newlevel, eventtime(tchg) selfloop name(levelchg)
restore
nwset
*# sp_change
nwdynam mail, inertia recip indeg alter(level) change(level levelchg)
*# sp_comp
* artificial membership, for illustration of present() and compchange() (not a recommended analysis):
* the employees who neither send nor receive an email before day 4.5 are absent at the start and join at day 4.5
quietly merge m:1 label using `early', keep(master match) nogenerate
gen byte member0 = seen == 1
tabulate member0
sort row
preserve
keep if member0 == 0
gen label2 = label
gen double tjoin = 4.5
gen byte joined = 1
nwset label label2 joined, eventtime(tjoin) selfloop name(membership)
restore
*# sp_compfit
nwdynam mail, inertia recip indeg present(member0) compchange(membership)
*# end
