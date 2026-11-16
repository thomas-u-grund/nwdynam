* Chapter 11: writes the event log that goldfish reads (ex/r/emails.csv with its actor list)
do setup.doh
nwwebuse manufacturing_nodes, nwclear
keep label
export delimited using "r/email_actors.csv", replace
nwwebuse manufacturing, nwclear
keep if t < 5
export delimited sender receiver t using "r/emails.csv", replace
*# end
