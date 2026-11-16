* Chapter 12: refusals and error messages
do setup.doh
mkmail 5
nwload mail, xvars
gen x = _n
*# li_rate
capture noisily nwdynam mail, submodel(rate) recip
display "return code: " _rc
*# li_choice_rate
capture noisily nwdynam mail, ego(x)
display "return code: " _rc
*# li_window
capture noisily nwdynam mail, inertiawindow(0)
display "return code: " _rc
*# li_intercept
capture noisily nwdynam mail, intercept
display "return code: " _rc
*# li_opp_rate
capture noisily nwdynam mail, submodel(rate) opportunities(x x)
display "return code: " _rc
*# li_net_rate
capture noisily nwdynam mail, submodel(rate) netterms(trans(mail))
display "return code: " _rc
*# li_coord_directed
capture noisily nwdynam mail, submodel(choice_coordination)
display "return code: " _rc
*# li_sub
capture noisily nwdynam mail, submodel(frequency)
display "return code: " _rc
*# li_noevent
nwset, mat((0,1\1,0)) name(plain) directed
capture noisily nwdynam plain
display "return code: " _rc
*# li_selfloop
clear
input sender receiver t
1 2 1
2 2 2
2 1 3
end
nwset sender receiver, directed eventtime(t) name(loops)
capture noisily nwdynam loops
display "return code: " _rc
*# li_restore
mkmail 5
quietly nwdynam mail, inertia recip indeg
estimates store m1
quietly nwdynam mail, inertia
estimates restore m1
capture noisily estat ic
display "return code: " _rc
*# li_two
quietly nwwebuse manufacturing_nodes, nwclear
rename label receiver
gen str5 division2 = "D" + division
keep receiver division2
tempfile rdiv
quietly save `rdiv'
quietly nwwebuse manufacturing, nwclear
keep if t < 5
quietly merge m:1 receiver using `rdiv', keep(match) nogenerate
keep sender division2 t
nwset sender division2, twomode eventtime(t) name(tm)
capture noisily nwdynam tm, trans
display "return code: " _rc
capture noisily nwdynam tm, inertiawindow(1)
display "return code: " _rc
*# end
