* Chapter 3 and 4: the first fit on the emails of a manufacturing company
do setup.doh
*# v_version
nwdynam version
*# first_data
nwwebuse manufacturing, nwclear
keep if t < 5
list in 1/5
count
display "span in days: " %5.2f t[_N] - t[1]
duplicates report sender receiver
*# first_isolates
quietly levelsof sender, local(a) clean
quietly levelsof receiver, local(b) clean
local seen : list a | b
display "employees in the log: " `: word count `seen''
tempfile mail
save `mail'
nwwebuse manufacturing_nodes, nwclear
quietly levelsof label, local(everyone) clean
display "employees in the company file: " `: word count `everyone''
local without : list everyone - seen
display "employees without an email: " `: word count `without''
local miss
foreach e of local without {
    local miss `miss',`e'
}
local miss = substr("`miss'", 2, .)
*# first_declare
use `mail', clear
nwset sender receiver, directed eventtime(t) name(mail)
nwaddnodes mail, nodenames(`miss')
nwset
*# first_choice
nwdynam mail
*# first_rate
nwdynam mail, submodel(rate)
*# first_replay
nwdynam, level(90)
*# first_ci
nwdynam, ci
*# first_ereturn
ereturn list
*# end
