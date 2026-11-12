cscript

* Event values as in goldfish (nwset sender receiver value, eventtime(time)): the value of an event is the increment of
* the tie (goldfish's increment column). A tie exists while its summed value is above 0 (unweighted effects); weighted
* effects use the summed value; a window keeps the summed value of the events inside the window. The reference numbers
* come from goldfish 1.6.12 (R 4.6) on the same 60 events (6 actors, two groups of simultaneous events):
*   vA: positive values 1..4;  vB: values with 0 and negative increments that dissolve ties (summed values never below 0).
* Every estimate must agree with goldfish to 1e-4 and the log likelihood to 1e-6 (relative).
*
* R code that made the events (written to ev.csv, listed below in the input block):
* set.seed(20261112)
* n <- 6; ne <- 60
* snd <- integer(ne); rcv <- integer(ne)
* for (k in 1:ne) { s <- sample(1:n, 1); r <- sample(setdiff(1:n, s), 1); snd[k] <- s; rcv[k] <- r }
* # bias some repetition
* for (k in seq(5, ne, by=4)) { snd[k] <- snd[k-3]; rcv[k] <- rcv[k-3] }
* tm <- cumsum(sample(1:3, ne, replace=TRUE))
* tm[c(10,11)] <- tm[9]; tm[c(31,32)] <- tm[30]   # two groups of simultaneous events
* tm <- sort(tm)
* # A: positive values 1..4
* vA <- sample(1:4, ne, replace=TRUE)
* # B: values with 0 and dissolving negatives, cumulative never below 0
* vB <- numeric(ne); cum <- matrix(0, n, n)
* for (k in 1:ne) {
*   c0 <- cum[snd[k], rcv[k]]
*   v <- if (c0 > 0 && k %% 4 == 0) -c0 else if (k %% 7 == 0) 0 else sample(1:2, 1)
*   vB[k] <- v; cum[snd[k], rcv[k]] <- c0 + v }
* x <- c(2.5, -1, 0.5, 3, 1, -2); g <- c(1, 2, 1, 2, 2, 1)
* d <- data.frame(snd, rcv, tm, vA, vB)
* write.csv(d, "ev.csv", row.names=FALSE)
* cat(paste(apply(d, 1, paste, collapse=" "), collapse="\n"), "\n")
*
* R code of the reference fits (prints coefficients, then the log likelihood):
* suppressMessages(library(goldfish))
* d <- read.csv("ev.csv"); n <- 6; lab <- paste0("A", 1:n)
* nodes <- defineNodes(data.frame(label=lab, present=TRUE, x=c(2.5,-1,0.5,3,1,-2), g=c(1,2,1,2,2,1), stringsAsFactors=FALSE))
* fit <- function(tag, form, sub, inc, directed=TRUE) {
*   ev <- data.frame(time=d$tm, sender=lab[d$snd], receiver=lab[d$rcv], increment=inc, stringsAsFactors=FALSE)
*   assign("ev", ev, globalenv()); assign("nodes", nodes, globalenv())
*   net <- defineNetwork(nodes=nodes, directed=directed); net <- linkEvents(x=net, changeEvent=ev, nodes=nodes)
*   assign("net", net, globalenv())
*   dep <- defineDependentEvents(events=ev, nodes=nodes, defaultNetwork=net); assign("dep", dep, globalenv())
*   m <- suppressWarnings(estimate(as.formula(paste("dep ~", form)), model="DyNAM", subModel=sub,
*          estimationInit=list(maxIterations=100, initialDamping=10)))
*   cat(tag, ": ", paste(sprintf("%.8f", c(coef(m), m$logLikelihood)), collapse=" "), "\n", sep="") }
* fit("A1", "inertia(net, weighted=TRUE)", "choice", d$vA)
* fit("A2", "inertia(net, weighted=TRUE) + recip(net, weighted=TRUE) + indeg(net, weighted=TRUE) + outdeg(net, weighted=TRUE)", "choice", d$vA)
* fit("A3", "indeg(net, weighted=TRUE) + outdeg(net, weighted=TRUE)", "rate", d$vA)
* fit("B1", "inertia(net) + recip(net) + indeg(net)", "choice", d$vB)
* fit("B2", "indeg(net) + outdeg(net)", "rate", d$vB)
* fit("B3", "inertia(net) + recip(net, weighted=TRUE) + trans(net) + alter(nodes$x)", "choice", d$vB)
* fit("C1", "inertia(net, weighted=TRUE, window=10) + recip(net, window=10) + same(nodes$g) + alter(nodes$x)", "choice", d$vA)
* fit("C2", "indeg(net, weighted=TRUE, window=8) + outdeg(net) + ego(nodes$x)", "rate", d$vA)
* fit("C3", "inertia(net, window=10) + indeg(net)", "choice", d$vB)
* fit("C4", "indeg(net, window=12) + outdeg(net, weighted=TRUE, window=12)", "rate", d$vB)
*
* Output (goldfish 1.6.12):
* A1: 0.05969556 -95.30813375
* A2: 0.04893543 -0.00503728 0.00138800 -0.04228915 -92.80953818
* A3: -0.00086567 -0.01083092 -107.36143380
* B1: 0.62774458 -0.12396545 0.12623936 -92.44713255
* B2: -0.04594828 0.20910414 -106.69708814
* B3: 0.89514471 -0.06614179 0.52476902 0.06689048 -90.46296811
* C1: 0.29634102 -0.15804754 0.30517480 0.06896162 -90.47074720
* C2: -0.13656306 -0.02429524 0.05479761 -105.69491420
* C3: 1.07970895 0.16300970 -90.03222704
* C4: -0.44739416 0.06900673 -103.73602452
* E1: fit("E1", "inertia(net) + recip(net) + indeg(net)", "choice", d$vB - 1)
*     0.41948555 0.09777224 0.00965213 -95.73599057
* E2: fit("E2", "indeg(net) + outdeg(net)", "rate", d$vB - 1)
*     -0.07475665 -0.00450931 -107.17592275
* and with directed=FALSE, subModel="choice_coordination":
* D1 (vB): fit("D1", "inertia(net) + indeg(net)", "choice_coordination", d$vB, directed=FALSE)
*     0.20874620 0.05003181 -161.44811731
* D0 (every value 1): 0.36279784 0.07180544 -160.24836654

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

* compare e(b) and e(ll) with goldfish (estimates in the order of e(b), then the log likelihood)
capture program drop chk
program chk
	args label
	local k = colsof(e(b))
	local ll : word `=`k'+2' of `0'
	forvalues j = 1/`k' {
		local ref : word `=`j'+1' of `0'
		if abs(el(e(b), 1, `j') - `ref') >= 1e-4 {
			di as error "`label': coefficient `j' " %12.8f el(e(b), 1, `j') " goldfish " %12.8f `ref'
			exit 9
		}
	}
	if reldif(e(ll), `ll') >= 1e-6 {
		di as error "`label': ll " %14.8f e(ll) " goldfish " %14.8f `ll'
		exit 9
	}
	di as text "  PASS: `label': ll " %12.6f e(ll) " (goldfish " %12.6f `ll' ")"
end

* actor covariates (by node label n1..n6)
capture program drop addcov
program addcov
	matrix __X = (2.5, -1, 0.5, 3, 1, -2)
	matrix __G = (1, 2, 1, 2, 2, 1)
	capture drop __id x g
	quietly generate __id = real(substr(_nwnode, 2, .))
	quietly generate x = __X[1, __id]
	quietly generate g = __G[1, __id]
	drop __id
end

clear
input snd rcv tm vA vB
2 1 1 2 2
5 2 3 3 2
5 4 4 1 1
5 2 5 3 -2
5 2 6 4 2
3 6 8 3 2
3 4 9 4 0
2 1 12 4 -2
3 6 14 1 2
1 4 14 3 1
5 6 14 4 2
4 3 19 3 2
1 4 21 4 2
2 6 24 2 0
1 4 25 1 1
5 2 26 4 -2
2 6 29 1 1
1 4 31 2 2
5 4 33 3 2
6 1 34 3 2
1 4 36 3 0
2 1 39 1 2
2 6 40 2 2
4 2 43 2 2
2 1 44 2 2
3 1 45 4 1
3 6 48 2 1
4 1 49 1 0
3 1 50 1 2
2 4 52 4 2
2 6 52 2 2
4 2 52 3 -2
2 4 59 4 2
4 6 60 3 1
3 6 63 4 0
2 1 65 4 -4
4 6 66 2 2
1 3 67 1 1
2 5 70 2 2
2 4 73 4 -4
1 3 76 3 2
1 5 78 3 0
3 5 80 3 1
1 6 83 1 1
1 5 86 4 1
6 4 89 2 2
5 1 91 2 2
4 2 93 2 2
6 4 95 3 0
6 3 97 2 2
4 6 98 2 1
3 1 99 1 -3
6 3 102 1 1
5 4 104 2 1
3 4 106 2 1
5 2 108 2 0
5 4 109 4 2
1 5 111 3 1
5 4 113 1 2
3 4 116 4 -1
end
tempfile ev
save `ev'

* ---- the event list carries the values (column 4); 1 without a value variable
nwclear
use `ev', clear
nwset snd rcv vA, eventtime(tm) name(a)
_nwsyntax a
mata: __evl = *(`netobj'->get_eventlist())
mata: assert(cols(__evl) == 4)
mata: assert(sum(__evl[., 4]) == 153)
nwclear
use `ev', clear
nwset snd rcv, eventtime(tm) name(u)
_nwsyntax u
mata: __evl = *(`netobj'->get_eventlist())
mata: assert(cols(__evl) == 4 & all(__evl[., 4] :== 1))
mata: mata drop __evl
di as text "  PASS: event list with values (column 4), 1 without a value variable"

* ---- without values = every value 1 (unchanged results)
nwdynam u
scalar ll_udef = e(ll)
nwdynam u, weightedinertia seed(1)
scalar ll_u = e(ll)
matrix b_u = e(b)
nwclear
use `ev', clear
generate one = 1
nwset snd rcv one, eventtime(tm) name(u1)
nwdynam u1, weightedinertia seed(1)
assert reldif(e(ll), ll_u) < 1e-12
assert mreldif(e(b), b_u) < 1e-12
di as text "  PASS: values all 1 = no value variable"

* ---- A: positive values, weighted effects
nwclear
use `ev', clear
nwset snd rcv vA, eventtime(tm) name(a)
nwdynam a, weightedinertia
chk "A1 choice weightedinertia" 0.05969556 -95.30813375
nwdynam a, weightedinertia weightedrecip weightedindeg weightedoutdeg
chk "A2 choice, all four weighted" 0.04893543 -0.00503728 0.00138800 -0.04228915 -92.80953818
nwdynam a, submodel(rate) weightedindeg weightedoutdeg
chk "A3 rate weighted indeg outdeg" -0.00086567 -0.01083092 -107.36143380
* positive values: the unweighted default fit is the same as without values
nwdynam a
assert reldif(e(ll), ll_udef) < 1e-12
di as text "  PASS: positive values leave the unweighted effects unchanged"

* ---- B: values 0 and below dissolve ties
nwclear
use `ev', clear
nwset snd rcv vB, eventtime(tm) name(b)
nwdynam b
chk "B1 choice default (native kernel)" 0.62774458 -0.12396545 0.12623936 -92.44713255
scalar ll_nat = e(ll)
nwdynam b, inertia recip indeg inertiawindow(1e9)
assert reldif(e(ll), ll_nat) < 1e-9
di as text "  PASS: B1 Mata evaluator = native evaluator with values"
nwdynam b, submodel(rate)
chk "B2 rate default (native kernel)" -0.04594828 0.20910414 -106.69708814
scalar ll_nat = e(ll)
nwdynam b, submodel(rate) indeg outdeg indegwindow(1e9)
assert reldif(e(ll), ll_nat) < 1e-9
di as text "  PASS: B2 Mata evaluator = native evaluator with values"
addcov
nwdynam b, inertia weightedrecip trans alter(x)
* e(b) order: inertia recip alter trans (goldfish printed trans before alter)
chk "B3 choice inertia weightedrecip trans alter" 0.89514471 -0.06614179 0.06689048 0.52476902 -90.46296811

* ---- E: tie values below 0 (vB - 1): inertia/recip need a value above 0, indeg/outdeg add the sign (-1)
nwclear
use `ev', clear
generate vC = vB - 1
nwset snd rcv vC, eventtime(tm) name(e)
nwdynam e
chk "E1 choice default, negative tie values (native kernel)" 0.41948555 0.09777224 0.00965213 -95.73599057
scalar ll_nat = e(ll)
nwdynam e, inertia recip indeg inertiawindow(1e9)
assert reldif(e(ll), ll_nat) < 1e-9
nwdynam e, submodel(rate)
chk "E2 rate default, negative tie values (native kernel)" -0.07475665 -0.00450931 -107.17592275
scalar ll_nat = e(ll)
nwdynam e, submodel(rate) indeg outdeg indegwindow(1e9)
assert reldif(e(ll), ll_nat) < 1e-9
di as text "  PASS: E Mata evaluator = native evaluator with negative tie values"

* ---- C: windows with values (and covariates)
nwclear
use `ev', clear
nwset snd rcv vA, eventtime(tm) name(c)
addcov
nwdynam c, weightedinertia inertiawindow(10) recipwindow(10) same(g) alter(x)
chk "C1 choice weighted window inertia, window recip, same, alter" 0.29634102 -0.15804754 0.30517480 0.06896162 -90.47074720
nwdynam c, submodel(rate) weightedindeg indegwindow(8) outdeg ego(x)
chk "C2 rate weighted window indeg, outdeg, ego" -0.13656306 -0.02429524 0.05479761 -105.69491420
nwclear
use `ev', clear
nwset snd rcv vB, eventtime(tm) name(c)
nwdynam c, inertiawindow(10) indeg
chk "C3 choice window inertia with dissolving values" 1.07970895 0.16300970 -90.03222704
nwdynam c, submodel(rate) indegwindow(12) weightedoutdeg outdegwindow(12)
chk "C4 rate window indeg, weighted window outdeg" -0.44739416 0.06900673 -103.73602452

* ---- D: undirected network, choice_coordination
nwclear
use `ev', clear
nwset snd rcv vB, eventtime(tm) undirected name(d)
nwdynam d, submodel(choice_coordination) inertia indeg
chk "D1 coordination with values" 0.20874620 0.05003181 -161.44811731
nwclear
use `ev', clear
nwset snd rcv, eventtime(tm) undirected name(d)
nwdynam d, submodel(choice_coordination) inertia indeg
chk "D0 coordination without values" 0.36279784 0.07180544 -160.24836654

* ---- a missing value is refused at nwset
nwclear
use `ev', clear
replace vA = . in 5
capture noisily nwset snd rcv vA, eventtime(tm) name(m)
assert _rc == 416
di as text "  PASS: missing event value refused (416)"

di as text "ALL PASS: test_nwdynam_eventvalues"
