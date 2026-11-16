* Chapter 8: opportunities(), oppframe() and event values
do setup.doh
* the events with their sequence number in the chronological order, and the actors with their division
quietly nwwebuse manufacturing, nwclear
keep if t < 5
gen long orig = _n
sort t orig
gen long evseq = _n
tempfile ev
quietly save `ev'
quietly nwwebuse manufacturing_nodes, nwclear
encode division, generate(div)
tempfile actors
quietly save `actors'
mkmail 5
*# op_ids
nwload mail, xvars
gen long actid = _n
gen label = _nwnode
quietly merge m:1 label using `actors', keep(master match) nogenerate
tempfile ids
quietly save `ids'
*# op_build
* the candidates of every event: the employees of the division of the sender, and the realised receiver
quietly use `ev', clear
rename sender label
quietly merge m:1 label using `ids', keep(match) nogenerate keepusing(div)
rename (label div) (sender sdiv)
rename sender label
tempfile evdiv
quietly save `evdiv'
quietly use `ids', clear
rename (label actid) (cand candid)
keep cand candid div
rename div sdiv
joinby sdiv using `evdiv'
keep evseq candid
tempfile opp1
quietly save `opp1'
quietly use `ev', clear
rename receiver label
quietly merge m:1 label using `ids', keep(match) nogenerate keepusing(actid)
rename actid candid
keep evseq candid
append using `opp1'
quietly duplicates drop
sort evseq candid
count
bysort evseq: gen n = _N
summarize n
keep evseq candid
tempfile opp
quietly save `opp'
*# op_fit_inmemory
quietly use `opp', clear
nwdynam mail, inertia recip indeg opportunities(evseq candid)
*# op_fit_frame
quietly use `ids', clear
frame create opps
frame opps: use `opp', clear
nwdynam mail, inertia recip indeg alter(level) opportunities(evseq candid) oppframe(opps)
*# op_value
* event values (artificial, to show the syntax): every email gets the value 1 or 2 at random
quietly nwwebuse manufacturing, nwclear
keep if t < 5
set seed 7
gen value = 1 + (runiform() < 0.3)
tabulate value
nwset sender receiver value, directed eventtime(t) name(valued)
nwdynam valued, inertia recip indeg
nwdynam valued, weightedinertia recip indeg
*# end
