* nwdynam 1.8.1: three defects found while writing the manual, with their goldfish 1.6.12 references (dev/dynam_fix181_crosscheck.R)
* on the first 631 emails of the manufacturing company (dev/fix181_*.csv):
*  1. netterms() with an event network (declared with eventtime()) changed the network's own event list: the second and every later fit that
*     used the same event network gave different, wrong estimates (and a sometimes unidentified term)
*  2. estat changepoints, plot stopped with "not found" and r(100) (the matrix was used after return matrix)
*  3. choice_coordination with a statistic of large values (indeg) stopped with "flat region", r(111), although the model is identified
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
quietly do `"`nwcore'/unw_core.do"'
do unw_dynam.do
do nwdynam.ado
do nwdynam_estat.ado

* ---- the emails as the event network mail (every one of the 154 employees is an actor)
clear
quietly import delimited using "dev/fix181_emails.csv", clear varnames(1)
assert _N == 631
quietly levelsof sender, local(a) clean
quietly levelsof receiver, local(b) clean
local seen : list a | b
tempfile emails
quietly save `emails'
quietly import delimited using "dev/fix181_actors.csv", clear varnames(1) stringcols(1)
quietly levelsof label, local(everyone) clean
quietly use `emails', clear
quietly nwset sender receiver, directed eventtime(t) name(mail)
local without : list everyone - seen
local miss
foreach e of local without {
	local miss `miss',`e'
}
if "`miss'" != "" quietly nwaddnodes mail, nodenames(`=substr("`miss'", 2, .)')
* the reporting hierarchy as an event network, every tie created at time 0
quietly import delimited using "dev/fix181_hier.csv", clear varnames(1) stringcols(1 2)
assert _N == 153
quietly nwset sender receiver, directed eventtime(t) name(hier)

* ---- 1. event network as the second network: goldfish ll -2297.03947783289, tie 0.971986964341019, inertia 1.715472042433106,
*         recip 2.893126700008462, indeg 0.110878617701675; the same in every fit (before 1.8.1 only the first fit was right)
_nwsyntax hier, other(fx_) max(1)
mata: HE0 = *(`fx_netobj'->get_eventlist())
forvalues k = 1/4 {
	quietly nwdynam mail, inertia recip indeg netterms(tie(hier))
	assert reldif(e(ll), -2297.03947783289) < 1e-6
	matrix b = e(b)
	assert reldif(b[1, 1], 1.715472042433106) < 1e-5 & reldif(b[1, 2], 2.893126700008462) < 1e-5
	assert reldif(b[1, 3], 0.110878617701675) < 1e-5 & reldif(b[1, 4], 0.971986964341019) < 1e-5
}
quietly nwdynam mail, inertia recip indeg netterms(tie(hier)) seed(7)
assert reldif(e(ll), -2297.03947783289) < 1e-6
* another fit in between, then the same fit again; the event list of hier is unchanged
quietly nwdynam mail, inertia recip indeg netterms(tie(hier) trans(hier))
quietly nwdynam mail, inertia recip indeg netterms(tie(hier))
assert reldif(e(ll), -2297.03947783289) < 1e-6
_nwsyntax hier, other(fx_) max(1)
mata: assert(HE0 == *(`fx_netobj'->get_eventlist()))
di as result "PASS: event network as second network, repeated fits"

* ---- 2. estat changepoints, plot (and outliers, plot) after a fit
quietly nwdynam mail, inertia recip indeg
estat changepoints, window(20) plot name(fix181_cp)
assert r(n_cpts) < . & r(N) == 631
capture estat changepoints, moment(variance) method(binseg) window(20) plot name(fix181_cp2)
assert _rc == 0
if r(n_cpts) > 0 {
	matrix C = r(cpts)
	assert colsof(C) == r(n_cpts) + 1
}
estat outliers, window(50) plot name(fix181_ol)
graph drop fix181_cp fix181_cp2 fix181_ol
di as result "PASS: estat changepoints, plot"

* ---- 3. coordination on the emails as undirected contacts (the actors that wrote or received an email): goldfish
*         indeg 0.09895999563 (ll -5250.956); inertia 1.87626305219, indeg 0.04106160058 (ll -4602.25); inertia 1.968821873 (ll -4621.636)
quietly import delimited using "dev/fix181_emails.csv", clear varnames(1)
quietly nwset sender receiver, undirected eventtime(t) name(contact)
quietly nwdynam contact, submodel(choice_coordination) indeg
assert reldif(el(e(b), 1, 1), 0.09895999563) < 5e-4 & reldif(e(ll), -5250.956) < 5e-4
quietly nwdynam contact, submodel(choice_coordination) inertia indeg
assert reldif(el(e(b), 1, 1), 1.87626305219) < 5e-4 & reldif(el(e(b), 1, 2), 0.04106160058) < 5e-4 & reldif(e(ll), -4602.25) < 5e-4
quietly nwdynam contact, submodel(choice_coordination) inertia
assert reldif(el(e(b), 1, 1), 1.968821873) < 5e-4 & reldif(e(ll), -4621.636) < 5e-4
di as result "PASS: coordination with indeg"
di as result "ALL TESTS PASSED: test_nwdynam_fix181"
