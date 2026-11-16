* Chapters 4 and 10: coordination sub-model, two-mode networks, opportunities, tied times, event values
do setup.doh
*# sm_coord_data
quietly nwwebuse manufacturing, nwclear
keep if t < 5
nwset sender receiver, undirected eventtime(t) name(contact)
nwset
*# sm_coord
nwdynam contact, submodel(choice_coordination) inertia
*# sm_coord2
nwdynam contact, submodel(choice_coordination) inertia trans
*# sm_two_data
* two-mode events: an employee (mode 1) writes to a division (mode 2: the division of the receiver, labelled D plus the
* identifier of the head of the division, so that a division is not mistaken for an employee)
quietly nwwebuse manufacturing_nodes, nwclear
rename label receiver
gen str5 recdiv = "D" + division
keep receiver recdiv
tempfile rdiv
quietly save `rdiv'
quietly nwwebuse manufacturing, nwclear
keep if t < 5
quietly merge m:1 receiver using `rdiv', keep(match) nogenerate
keep sender recdiv t
rename recdiv division
nwset sender division, twomode eventtime(t) name(tm)
nwset
*# sm_two
nwdynam tm
nwdynam tm, submodel(rate) outdeg
*# sm_ties
mkmail 5
nwdynam mail, inertia recip indeg
quietly nwwebuse manufacturing, nwclear
keep if t < 5
replace t = floor(t * 4) / 4
nwset sender receiver, directed eventtime(t) name(coarse)
nwdynam coarse, inertia recip indeg
nwdynam coarse, inertia recip indeg sequential
*# end
