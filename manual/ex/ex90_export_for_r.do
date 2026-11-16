* Chapter 11: writes the two event logs that goldfish reads (ex/r/emails.csv with its actor list, ex/r/socevo.csv)
do setup.doh
nwwebuse manufacturing_nodes, nwclear
keep label
export delimited using "r/email_actors.csv", replace
nwwebuse manufacturing, nwclear
keep if t < 5
export delimited sender receiver t using "r/emails.csv", replace
import delimited using "data/socialevolution_calls.csv", clear
rename eventtime t
export delimited sender receiver t using "r/socevo.csv", replace
*# end
