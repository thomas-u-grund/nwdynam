cscript

* Long messages wrap at word boundaries (1.2.3). Stata cuts a display line that is longer than the linesize at exactly that
* column, in the middle of a word; the messages are now SMCL paragraphs. The test sets linesize 80 (as in the book),
* triggers two errors and the note on unidentified effects, and checks the log text.

local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
adopath ++ `"`nwcore'"'
capture adopath ++ `"`nwcore'/lib"'
adopath ++ `"`c(pwd)'"'
adopath ++ `"`c(pwd)'/lib"'
adopath ++ `"`c(pwd)'/lib/plugins/macos"'

* _mw_check logfile "message": the message must appear in the log, joined from its output lines (trimmed, joined with one
* space) exactly as written, on at least two lines, none longer than the linesize and none a hard-wrap continuation ("> ").
capture program drop _mw_check
program _mw_check
	args logf msg
	local start = substr(`"`macval(msg)'"', 1, 30)
	tempname fh
	local found 0
	local nl 0
	local hard 0
	local toolong 0
	local joined ""
	file open `fh' using `"`logf'"', read text
	file read `fh' line
	while r(eof) == 0 {
		if `found' == 1 & strlen(`"`macval(joined)'"') < strlen(`"`macval(msg)'"') {
			if substr(`"`macval(line)'"', 1, 2) == "> " local hard 1
			if strlen(`"`macval(line)'"') > 80 local toolong 1
			local t = trim(`"`macval(line)'"')
			local joined `"`macval(joined)' `macval(t)'"'
			local ++nl
		}
		if `found' == 0 & strpos(trim(`"`macval(line)'"'), `"`macval(start)'"') == 1 {
			local found 1
			if strlen(`"`macval(line)'"') > 80 local toolong 1
			local joined = trim(`"`macval(line)'"')
			local nl 1
		}
		file read `fh' line
	}
	file close `fh'
	di as text "  lines: `nl'  hard wrap: `hard'  too long: `toolong'"
	assert `found' == 1
	assert `nl' >= 2
	assert `hard' == 0
	assert `toolong' == 0
	assert `"`macval(joined)'"' == `"`macval(msg)'"'
	di as text "  PASS: wraps at word boundaries: " substr(`"`macval(msg)'"', 1, 50) "..."
end

set linesize 80
tempfile l1 l2 l3

* 1. error: a network that is not an event network
nwclear
clear
mata: rseed(4); R = runiform(10, 10) :< 0.2; R = R :* (1 :- I(10)); st_matrix("A", R)
nwset, mat(A) name(static) nooutput
log using `l1', text replace name(mw1)
capture noisily nwdynam static
local rc = _rc
log close mw1
assert `rc' == 198
_mw_check `l1' "nwdynam requires a network declared via nwset's eventtime() option (an event-type temporal network) - see nwset."

* 2. error: a window that is not positive; 3. note: effects that are not identified (toy log of ten events)
nwclear
clear
input sender receiver t
1 2 1
1 3 2
2 1 3
1 2 4
3 2 5
2 3 6
1 3 7
3 1 8
2 1 9
1 3 10
end
nwset sender receiver, eventtime(t) directed name(chat)
log using `l2', text replace name(mw2)
capture noisily nwdynam chat, submodel(choice) inertia inertiawindow(0)
local rc = _rc
log close mw2
assert `rc' == 198
_mw_check `l2' "inertiawindow() must be a positive number (in the same time units as the network's own eventtime())."

log using `l3', text replace name(mw3)
nwdynam chat
log close mw3
assert `"`e(unidentified)'"' != ""
_mw_check `l3' "These effects are probably not identified (too few events, a perfectly predicted pattern, or effects that are linear functions of each other). Do not interpret their estimates; drop or recode the effect."

di as result "ALL TESTS PASSED"
