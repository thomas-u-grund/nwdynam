* Default calls of the plots of the released nwdynam 1.9.1 (estat outliers, plot / estat changepoints, plot); each default picture is
* exported to the directory in the global NWP_OUT. Run with the released package and with the changed one: the PNG files must be identical.
* Used by test_nwdynam_plots.do (reference directory: environment variable NWDYNAM_PLOT_REF) and by hand. cscript first: it fixes the rendering state.
set graphics off
nwclear
capture mkdir "$NWP_OUT"
set seed 20261120
clear
set obs 60
generate long sender = ceil(6 * runiform())
generate long receiver = ceil(5 * runiform())
replace receiver = receiver + (receiver >= sender)
generate double eventtime = _n + runiform() / 4
sort eventtime
nwset sender receiver, eventtime(eventtime) name(toy) directed dropvars
nwdynam toy, submodel(choice) inertia
estat outliers, method(iqr) parameter(1) plot
quietly graph export "$NWP_OUT/outliers.png", replace width(600)
graph drop _all
estat changepoints, window(3) plot
quietly graph export "$NWP_OUT/changepoints.png", replace width(600)
graph drop _all
estat outliers, method(iqr) parameter(1) plot name(gn1)
quietly graph export "$NWP_OUT/outliers_named.png", name(gn1) replace width(600)
graph drop _all
