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
estat outliers, window(50) plot name(ol) export(../fig/outliers.pdf) xsize(5.5) ysize(3.4)
*# po_choiceplot
estat choiceplot, topk(1 3 5 10) name(cp) export(../fig/choiceplot.pdf) xsize(5.5) ysize(3.6) frame(cpf)
return list
*# po_rateint
nwdynam mail, submodel(rate) intercept indeg outdeg
estat eventlik, frame(evr)
return list
*# po_rateplot
nwdynam mail, submodel(rate)
estat rateplot, name(ra) export(../fig/rateplot_actor.pdf) xsize(4.5) ysize(4.4) frame(raf)
return list
nwdynam mail, submodel(rate) intercept indeg outdeg
estat rateplot, over(time) periods(5) topactors(4) name(rt) export(../fig/rateplot_time.pdf) xsize(5.5) ysize(3.6)
*# end
