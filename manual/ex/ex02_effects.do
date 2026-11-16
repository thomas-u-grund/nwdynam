* Chapter 10: effects, covariates, comparing models
do setup.doh
*# eff_attr
quietly nwwebuse manufacturing_nodes, nwclear
encode division, generate(div)
tempfile actors
quietly save `actors'
mkmail 5
nwload mail, xvars
gen long row = _n
gen label = _nwnode
merge m:1 label using `actors', keep(master match) nogenerate
sort row
list label level div in 1/4
tabulate level
*# eff_base
nwdynam mail, inertia recip indeg
estimates store base
*# eff_same
nwdynam mail, inertia recip indeg same(div)
estimates store samediv
estimates table base samediv, b(%9.3f) stats(ll N)
lrtest samediv base
*# eff_closure
nwdynam mail, inertia recip indeg trans
nwdynam mail, inertia recip indeg nodetrans
*# eff_window
nwdynam mail, inertiawindow(0.5) recip indeg
nwdynam mail, weightedinertia recip indeg
*# eff_rate
nwdynam mail, submodel(rate) indeg outdeg
estimates store rate0
nwdynam mail, submodel(rate) indeg outdeg ego(level)
estimates store rate1
lrtest rate1 rate0
*# eff_ic
estat ic
*# end
