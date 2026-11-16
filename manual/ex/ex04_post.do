* Chapter 9: postestimation
do setup.doh
mkmail 5
nwdynam mail, inertia recip indeg
*# po_ic
estat ic
return list
*# po_eventlik
estat eventlik, frame(evl)
frame evl: list in 1/5
frame evl: summarize ll
*# po_eventprob
estat eventprob, frame(evp)
frame evp: list in 1/6
frame evp: summarize p
*# po_outliers
estat outliers, method(top) parameter(5) frame(out)
return list
*# po_outliers2
estat outliers, window(50)
matrix list r(outliers)
*# po_cpts
estat changepoints, window(20)
return list
*# po_cpts2
estat changepoints, moment(variance) method(binseg) window(20)
matrix list r(cpts)
*# po_plot
quietly estat outliers, window(50) plot name(ol)
graph export "../fig/outliers.pdf", name(ol) replace
*# po_rateint
nwdynam mail, submodel(rate) intercept indeg outdeg
estat eventlik, frame(evr)
return list
*# end
