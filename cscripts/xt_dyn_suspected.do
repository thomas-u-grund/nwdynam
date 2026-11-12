cscript

* xt_dyn_suspected.do -- DRAFT, not yet run.
* Checks of SUSPECTED bugs and code/help or nwdynam/goldfish mismatches found by reading the code (details and
* file:line in xt_NOTES.md). Each check runs on its own (capture noisily { }); the file reports every failing check
* by name and stops with r(9) at the end if any failed. A failure here is a finding to look at, not necessarily a
* regression. Checks marked DESIGN encode a behaviour that is arguably right but not promised by the help file.
* The reference fit of the invariance checks is the W40 default choice model (goldfish: 1.0497675 1.0747602
* -0.2707421, ll -38.1300248; xt_dyn_edges.do).
* Expected runtime: about 30 seconds (Stata 17 BE).
* Run from the repository root; the test runner sets the adopath (core, this repository, the plugin folder).
* No sysdir/adopath here.

capture program drop nwinstall
local nwcore : env NWCORE
if `"`nwcore'"' != "" quietly do `"`nwcore'/unw_core.do"'
quietly do unw_dynam.do

mata: XT_W40 = (1,3 \ 1,3 \ 1,3 \ 1,3 \ 2,1 \ 4,2 \ 2,4 \ 2,4 \ 1,3 \ 4,3 \ 3,1 \ 4,2 \ 3,2 \ 4,2 \ 1,3 \ 1,3 \ 2,4 \ 4,2 \ 2,4 \ 2,4 \ 3,1 \ 3,1 \ 1,2 \ 1,3 \ 3,4 \ 4,1 \ 1,3 \ 3,1 \ 1,2 \ 1,3 \ 4,2 \ 2,1 \ 4,2 \ 4,2 \ 3,2 \ 2,3 \ 1,3 \ 2,4 \ 4,2 \ 2,3)

capture program drop xt_w40
program xt_w40
	clear
	quietly set obs 40
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store(., ("sender", "receiver"), XT_W40)
	quietly gen double t = _n
end

local fails 0
local failed ""

nwclear
xt_w40
nwset sender receiver, eventtime(t) directed name(w40)
nwdynam w40, seed(1)
assert reldif(e(ll), -38.1300248) < 1e-5
matrix XB = e(b)
scalar XLL = e(ll)

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CORE BUG S1 (cause in core nwset; nwdynam refuses the events with r(416) until the core is fixed; ids of 8+ digits already fail: %9.0g prints 1.00e+07)
di as text "S1: numeric actor ids above 16,777,216 (not exact in float)"
di as text "{hline 60}"
* core nwset.ado:285-286 labels events "n" + strofreal(id) (%9.0g); nwfromedge.ado:218, 246 builds the node names
* from a float dictionary with %18.0g: labels differ, events get actor index 0.
capture noisily {
	nwclear
	xt_w40
	replace sender = 20000000 + 2*sender - 1
	replace receiver = 20000000 + 2*receiver - 1
	nwset sender receiver, eventtime(t) directed name(wbig)
	nwdynam wbig, seed(1)
	assert e(nodes) == 4 & e(N) == 40
	assert mreldif(e(b), XB) < 1e-10
}
if _rc == 416 di as text "  CORE BUG S1 (not nwdynam; nwdynam refuses the unmatched events cleanly, r(416))"
else if _rc {
	local ++fails
	local failed "`failed' S1"
	di as error "  FAIL S1 (r(`=_rc')): large numeric ids"
}
else di as text "  PASS: S1"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CORE BUG S2 (same cause as S1)
di as text "S2: numeric actor ids 1e9 .. 4e9"
di as text "{hline 60}"
capture noisily {
	nwclear
	xt_w40
	recast double sender receiver
	replace sender = sender * 1e9
	replace receiver = receiver * 1e9
	nwset sender receiver, eventtime(t) directed name(wbil)
	nwdynam wbil, seed(1)
	assert e(nodes) == 4 & e(N) == 40
	assert mreldif(e(b), XB) < 1e-10
}
if _rc == 416 di as text "  CORE BUG S2 (not nwdynam; nwdynam refuses the unmatched events cleanly, r(416))"
else if _rc {
	local ++fails
	local failed "`failed' S2"
	di as error "  FAIL S2 (r(`=_rc')): ids of 1e9 and more"
}
else di as text "  PASS: S2"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S3 (native choice/rate engines fit the event with sender index 0 silently, e(N)=41; multi engine r(3301); missing receiver: silent on W40, r(111) on small data)
di as text "S3: an event row with a missing sender"
di as text "{hline 60}"
* the event gets actor index 0 (label "n." not found); DynamState indexes column 0: expected a raw Mata error
capture noisily {
	nwclear
	xt_w40
	set obs 41
	replace receiver = 2 in 41
	replace t = 41 in 41
	* the core now refuses the event at nwset (r(416), a message, nothing created)
	capture noisily nwset sender receiver, eventtime(t) directed name(wmis)
	assert _rc == 416
}
if _rc {
	local ++fails
	local failed "`failed' S3"
	di as error "  FAIL S3: missing sender gives a raw Mata error or is counted as an event"
}
else di as text "  PASS: S3"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S4 (choice: timeless event used as last event; rate+intercept: r(111) after "initial values not feasible")
di as text "S4: an event row with a missing time"
di as text "{hline 60}"
* nwset stores the missing time; sorting puts it last: the event silently becomes the LAST event (and with
* intercept the waiting time is missing)
capture noisily {
	nwclear
	xt_w40
	replace t = . in 10
	* the core now refuses the event at nwset (r(416), a message, nothing created)
	capture noisily nwset sender receiver, eventtime(t) directed name(wmt)
	assert _rc == 416
}
if _rc {
	local ++fails
	local failed "`failed' S4"
	di as error "  FAIL S4: an event without a time is used (placed last)"
}
else di as text "  PASS: S4"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S5 (accepted, ll = -1e300; goldfish refuses)
di as text "S5: opportunities() that do not contain the realised receiver"
di as text "{hline 60}"
* goldfish 1.6.12 refuses ("Active node 3 not available in event 1"); nwdynam.ado:686-735 checks only the ranges
capture noisily {
	nwclear
	xt_w40
	nwset sender receiver, eventtime(t) directed name(w40)
	clear
	quietly set obs 160
	gen long oev = ceil(_n / 4)
	gen long oact = _n - 4 * (oev - 1)
	* event 1 is 1 -> 3: remove actor 3 from its opportunities
	drop if oev == 1 & oact == 3
	capture noisily nwdynam w40, inertia opportunities(oev oact) seed(1)
	assert _rc != 0
	assert _rc < 3000 | _rc > 3999
}
if _rc {
	local ++fails
	local failed "`failed' S5"
	di as error "  FAIL S5: realised receiver outside the opportunity set accepted (goldfish refuses)"
}
else di as text "  PASS: S5"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S6 (raw Mata r(3200) conformability error)
di as text "S6: opportunities() data with missing values"
di as text "{hline 60}"
* the range checks skip missing values (nwdynam.ado:694, 699 'if !missing()'); the Mata matrix then holds a missing
* subscript. Wanted: the rows ignored (= fit without them) or a clear refusal.
capture noisily {
	clear
	quietly set obs 160
	gen long oev = ceil(_n / 4)
	gen long oact = _n - 4 * (oev - 1)
	quietly set obs 162
	replace oact = 2 in 161
	replace oev = 7 in 162
	capture noisily nwdynam w40, inertia opportunities(oev oact) seed(1)
	local rc = _rc
	if `rc' == 0 assert reldif(e(ll), -41.213570853) < 1e-5
	else assert `rc' < 3000 | `rc' > 3999
}
if _rc {
	local ++fails
	local failed "`failed' S6"
	di as error "  FAIL S6: missing values in the opportunities() data give a raw Mata error"
}
else di as text "  PASS: S6"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S7 (cosmetic: scalar leaks)
di as text "S7: scalar __nwdynam_sequential left behind after an error"
di as text "{hline 60}"
* nwdynam.ado:656-657 creates the scalar before the event-count check (659-662) and the opportunities()/tie()/two-mode
* checks (686-772); only the success path drops it (816).
capture noisily {
	nwclear
	clear
	input long sender long receiver double t
	1 2 1
	end
	nwset sender receiver, eventtime(t) directed name(one)
	capture scalar drop __nwdynam_sequential
	capture nwdynam one
	assert _rc == 2001
	capture confirm scalar __nwdynam_sequential
	assert _rc != 0
}
if _rc {
	local ++fails
	local failed "`failed' S7"
	di as error "  FAIL S7: scalar __nwdynam_sequential remains after an error"
}
else di as text "  PASS: S7"
capture frame change default
capture frame drop __nwdynam_native
capture scalar drop __nwdynam_sequential

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S8 (cosmetic: V = 0 not flagged; e(unidentified) empty)
di as text "S8: an effect that is constant zero is not reported as unidentified"
di as text "{hline 60}"
* inertiawindow(0.5) on integer times: no earlier event is ever within 0.5, the statistic is 0 for every candidate,
* the likelihood is flat: ll = null = -40 ln 3 = -43.9444915467. The check nwdynam.ado:855-860 flags only a
* missing or > 100 standard error; nwrem also flags a variance <= 0.
capture noisily {
	nwclear
	xt_w40
	nwset sender receiver, eventtime(t) directed name(w40)
	nwdynam w40, inertiawindow(0.5) seed(1)
	assert abs(e(ll) - (-43.9444915467)) < 1e-6
	assert strpos(" `e(unidentified)' ", " inertia ") > 0
}
if _rc {
	local ++fails
	local failed "`failed' S8"
	di as error "  FAIL S8: flat effect not flagged in e(unidentified)"
}
else di as text "  PASS: S8"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S9 (minor: every window option (-1) silently means no window)
di as text "S9: inertiawindow(-1) is the internal 'not given' value"
di as text "{hline 60}"
* nwdynam.ado:16 INERTIAWINDOW(real -1); 362: -1 = not given. inertiawindow(-2) is refused (r(198)), inertiawindow(-1)
* silently fits the default model (inertia recip indeg).
capture noisily {
	capture noisily nwdynam w40, inertiawindow(-1) seed(1)
	assert _rc == 198
}
if _rc {
	local ++fails
	local failed "`failed' S9"
	di as error "  FAIL S9: inertiawindow(-1) accepted as 'no window'"
}
else di as text "  PASS: S9"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S10 (minor: seed(-1) silently unseeded)
di as text "S10: seed(-1) is the internal 'not given' value"
di as text "{hline 60}"
* nwdynam.ado:16 SEED(integer -1), 642-644: -1 means 'no seed'. The 8-effect set runs on the multi engine with random
* starts; with seed(-1) the result depends on the current random-number state.
capture noisily {
	set seed 5
	capture noisily nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(-1)
	local rc1 = _rc
	if `rc1' == 0 {
		matrix XS1 = e(b)
		set seed 6
		nwdynam w40, inertia recip indeg outdeg trans cycle commonsender commonreceiver seed(-1)
		assert mreldif(e(b), XS1) == 0
	}
}
if _rc {
	local ++fails
	local failed "`failed' S10"
	di as error "  FAIL S10: seed(-1) is ignored"
}
else di as text "  PASS: S10"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S11 (alter()/ego(): r(111) after "initial values not feasible"; same(): missing silently treated as its own category)
di as text "S11: covariate with a missing value"
di as text "{hline 60}"
* same()/diff()/sim()/alter()/ego()/tertius()/egoalterint() are read with st_data (nwdynam.ado:603-640) and never
* checked for missing values.
capture noisily {
	nwload w40, xvars
	gen double x = _n
	replace x = . in 2
	// wanted: a clean refusal (198 or 416), not r(111) after optimizer messages; same() must not run silently
	capture noisily nwdynam w40, alter(x) seed(1)
	assert inlist(_rc, 198, 416)
	capture noisily nwdynam w40, submodel(rate) ego(x) seed(1)
	assert inlist(_rc, 198, 416)
	capture noisily nwdynam w40, same(x) seed(1)
	assert inlist(_rc, 198, 416)
}
if _rc {
	local ++fails
	local failed "`failed' S11"
	di as error "  FAIL S11: missing covariate value accepted"
}
else di as text "  PASS: S11"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* S12 (FIXED): goldfish adds the event values as increments of weighted statistics (R 2026-11-12:
* 0.1389639466 / -39.2302174436; rate -0.009244684541 -0.038039028077 / -55.218661708709). The core event list now
* carries the value of each event and nwdynam adds it (counting events gave 0.2889076707 / -37.7174439566).
di as text "S12: weighted effects on a valued event network use the event values (goldfish)"
di as text "{hline 60}"
capture noisily {
	nwclear
	xt_w40
	gen double w = mod(_n, 3) + 1
	nwset sender receiver w, eventtime(t) directed name(wval)
	nwdynam wval, weightedinertia seed(1)
	assert abs(el(e(b),1,1) - 0.1389639466) < 1e-4
	assert reldif(e(ll), -39.2302174436) < 1e-6
	nwdynam wval, submodel(rate) weightedindeg weightedoutdeg seed(1)
	assert abs(el(e(b),1,1) - (-0.009244684541)) < 1e-4
	assert abs(el(e(b),1,2) - (-0.038039028077)) < 1e-4
	assert reldif(e(ll), -55.218661708709) < 1e-6
}
if _rc {
	local ++fails
	local failed "`failed' S12"
	di as error "  FAIL S12: weighted effects on a valued event network"
}
else di as text "  PASS: S12 (event values added as in goldfish)"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
* CONFIRMED BUG S13 (owner 2026-11-12: covariates must be matched via _nwnode when present, not by row)
di as text "S13 (DESIGN): covariate dataset re-sorted although _nwnode is there"
di as text "{hline 60}"
* read by row position (help: 'one row per actor, in actor order'). Closed form for alter(x4), x4 = (0,0,0,1):
* -0.496436886314, ll -43.2266914797 (xt_dyn_closedform.do).
capture noisily {
	nwclear
	xt_w40
	nwset sender receiver, eventtime(t) directed name(w40)
	nwload w40, xvars
	gen double x4 = (_n == 4)
	gsort -_nwnode
	capture noisily nwdynam w40, alter(x4) seed(1)
	if _rc == 0 assert abs(el(e(b),1,1) - (-0.496436886314)) < 1e-4
}
if _rc {
	local ++fails
	local failed "`failed' S13"
	di as error "  FAIL S13: covariate read by row position from a re-sorted dataset"
}
else di as text "  PASS: S13"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
* CONFIRMED BUG N2 (triage): self-loop refusal (unw_dynam.do:167) passes a ~140-character text to Mata _error(), which
* fails with r(3300) "argument out of range"; the message is lost. Wanted r(3498) with the message.
di as text "{hline 60}"
di as text "N2: self-loop event gives r(3498) with its message"
di as text "{hline 60}"
capture noisily {
	nwclear
	xt_w40
	replace receiver = sender in 3
	nwset sender receiver, eventtime(t) directed name(wloop)
	capture noisily nwdynam wloop, seed(1)
	assert _rc == 3498
}
if _rc {
	local ++fails
	local failed "`failed' N2"
	di as error "  FAIL N2: self-loop refusal loses its message (r(3300))"
}
else di as text "  PASS: N2"
capture frame change default
capture frame drop __nwdynam_native

* CONFIRMED BUG I11 (notes item 11, triage): when optimize() fails inside the native unit1 fit, the user is left in
* frame __nwdynam_native (own data apparently gone) and every later nwdynam call fails with r(3598) "frame
* __nwdynam_native already defined" until the frame is dropped by hand (unw_dynam.do:457-460, 519-520, 683-684).
* Fault injection: 9 events, the last one with a missing receiver: the default choice fit stops with r(111).
di as text "{hline 60}"
di as text "I11: frame __nwdynam_native after a failed native fit"
di as text "{hline 60}"
capture noisily {
	nwclear
	clear
	quietly set obs 9
	quietly gen long sender = .
	quietly gen long receiver = .
	mata: st_store((1::8), ("sender","receiver"), (1,2\2,3\3,1\1,3\2,1\4,2\1,4\3,4))
	quietly replace sender = 2 in 9
	quietly gen double t = _n
	* the core now refuses the event with the missing receiver at nwset (r(416)), so this fault no longer reaches the
	* native fit; nwdynam's own checks also come before the helper frame is created
	capture noisily nwset sender receiver, eventtime(t) directed name(wfail)
	assert _rc == 416
	assert "`c(frame)'" == "default"
	nwclear
	xt_w40
	nwset sender receiver, eventtime(t) directed name(w40)
	capture noisily nwdynam w40, seed(1)
	assert _rc == 0
}
if _rc {
	local ++fails
	local failed "`failed' I11"
	di as error "  FAIL I11: failed native fit leaves frame __nwdynam_native (now in frame `c(frame)')"
}
else di as text "  PASS: I11"
capture frame change default
capture frame drop __nwdynam_native

* ---------------------------------------------------------------------------------------------------------------------
di as text "{hline 60}"
di as text "N1: replay (nwdynam alone, nwdynam, level(90), estimates replay)"
di as text "{hline 60}"
capture noisily {
	nwclear
	xt_w40
	nwset sender receiver, eventtime(t) directed name(wrp)
	nwdynam wrp, seed(1)
	matrix XR = e(b)
	estimates store xtrp
	nwdynam
	assert mreldif(e(b), XR) == 0
	nwdynam, level(90)
	ereturn clear
	estimates replay xtrp
	estimates restore xtrp
	nwdynam
	assert mreldif(e(b), XR) == 0
	estimates drop xtrp
}
if _rc {
	local ++fails
	local failed "`failed' N1"
	di as error "  FAIL N1: replay"
}
else di as text "  PASS: N1"

if `fails' > 0 {
	di as error "suspected problems confirmed: `failed'"
	exit 9
}
di "ALL PASS"
