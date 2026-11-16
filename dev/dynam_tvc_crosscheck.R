# (terms in the order of nwdynam's effect list, so that e(b) and goldfish's coefficients line up)
# nwdynam time-varying covariates against goldfish (R 4.x, goldfish 1.6.12).
# Usage: Rscript dev/dynam_tvc_crosscheck.R <output directory>
# For every case it writes events_<case>.csv (time, sender, receiver), nodes_<case>.csv (label, x, y, g: the values at the start),
# changes_<var>_<case>.csv (time, node, value; value = the new value, or the increment for the increment cases) and one file
# results.csv (case, term, estimate, se, ll) with the goldfish fit. The Stata side (dev/dynam_tvc_crosscheck.do) reads them.
suppressMessages(library(goldfish))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# one case: n actors, N events, change events of the variables in `chg` (character vector of x, y, g), `inc` = increments
mkcase <- function(seed, n, N, chg, inc = FALSE, directed = TRUE, nchg = 15, ties = TRUE, edge = TRUE) {
  set.seed(seed)
  lab <- paste0("A", 1:n)
  snd <- integer(N); rcv <- integer(N)
  for (k in 1:N) {
    s <- sample(1:n, 1); r <- sample(setdiff(1:n, s), 1); snd[k] <- s; rcv[k] <- r
    if (!directed) { }
  }
  # some repetition
  for (k in seq(5, N, by = 4)) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }
  tm <- round(cumsum(runif(N, 0.3, 3)), 2)
  if (ties) { for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1] }
  tm <- sort(tm)
  ev <- data.frame(time = tm, sender = lab[snd], receiver = lab[rcv], increment = 1, stringsAsFactors = FALSE)
  nodes <- data.frame(label = lab, present = TRUE, x = round(rnorm(n), 2), y = round(runif(n, 0, 2), 2),
                      g = sample(1:3, n, replace = TRUE), stringsAsFactors = FALSE)
  chl <- list()
  for (v in chg) {
    ct <- sort(round(runif(nchg, min(tm) - 2, max(tm) + 2), 2))
    # changes exactly at event times (ties between a change and a dependent event) and two at the same time for one actor
    ct[2:4] <- tm[c(3, 10, 20)]
    ct[5:6] <- tm[15]
    ct <- sort(ct)
    cn <- sample(lab, nchg, replace = TRUE)
    if (v == "g") cv <- sample(1:3, nchg, replace = TRUE) else cv <- round(rnorm(nchg), 2)
    if (inc && v != "g") cv <- round(rnorm(nchg, 0, 0.5), 2)
    if (!edge) { ct <- sort(round(runif(nchg, min(tm) + 0.01, max(tm) - 0.01), 2)); ct[2:3] <- tm[c(3, 10)] ; ct <- sort(ct) }
    chl[[v]] <- data.frame(time = ct, node = cn, value = cv, stringsAsFactors = FALSE)
  }
  list(lab = lab, ev = ev, nodes = nodes, chl = chl, inc = inc, directed = directed)
}

fitgf <- function(cs, form, sub, useChanges = TRUE) {
  nodes <- defineNodes(cs$nodes)
  assign("ev", cs$ev, globalenv())
  if (useChanges) for (v in names(cs$chl)) {
    d <- cs$chl[[v]]
    if (cs$inc) names(d)[3] <- "increment" else names(d)[3] <- "replace"
    nm <- paste0("chg_", v)
    assign(nm, d, globalenv())
    nodes <- eval(substitute(linkEvents(nodes, changeEvents = CH, attribute = V), list(CH = as.name(nm), V = v)))
  }
  assign("nodes", nodes, globalenv())
  net <- defineNetwork(nodes = nodes, directed = cs$directed)
  net <- linkEvents(net, changeEvent = ev, nodes = nodes)
  assign("net", net, globalenv())
  dep <- defineDependentEvents(events = ev, nodes = nodes, defaultNetwork = net)
  assign("dep", dep, globalenv())
  m <- suppressWarnings(estimate(as.formula(paste("dep ~", form)), model = "DyNAM", subModel = sub,
          estimationInit = list(maxIterations = 200, initialDamping = 10, maxScoreStopCriterion = 1e-9)))
  list(b = coef(m), se = m$standardErrors, ll = as.numeric(m$logLikelihood))
}

cases <- list(
  list(id = "c1", seed = 11, n = 6, N = 40, chg = "x", form = "inertia(net) + alter(nodes$x)", sub = "choice"),
  list(id = "c2", seed = 12, n = 7, N = 60, chg = c("x", "g"), form = "inertia(net) + same(nodes$g) + diff(nodes$x)", sub = "choice"),
  list(id = "c3", seed = 13, n = 7, N = 60, chg = c("x", "y"), form = "inertia(net) + alter(nodes$y) + egoAlterInt(nodes$x, nodes$y)", sub = "choice"),
  list(id = "c4", seed = 14, n = 7, N = 70, chg = "x", form = "inertia(net) + indeg(net) + tertius(net, nodes$x)", sub = "choice"),
  list(id = "c5", seed = 15, n = 7, N = 60, chg = "x", form = "indeg(net) + ego(nodes$x)", sub = "rate"),
  list(id = "c6", seed = 16, n = 7, N = 70, chg = "x", form = "outdeg(net) + tertius(net, nodes$x)", sub = "rate"),
  list(id = "c7", seed = 17, n = 7, N = 60, chg = "x", form = "1 + indeg(net) + ego(nodes$x)", sub = "rate"),
  list(id = "c8", seed = 18, n = 7, N = 50, chg = c("x", "g"), form = "inertia(net) + same(nodes$g) + sim(nodes$x) + alter(nodes$x) + tertius(net, nodes$x)", sub = "choice_coordination", directed = FALSE),
  list(id = "c9", seed = 19, n = 7, N = 60, chg = c("x", "y"), inc = TRUE, form = "inertia(net) + diff(nodes$x) + alter(nodes$y) + egoAlterInt(nodes$x, nodes$y)", sub = "choice"),
  list(id = "c10", seed = 20, n = 7, N = 70, chg = "x", form = "inertia(net, window = 8) + recip(net) + sim(nodes$x) + alter(nodes$x)", sub = "choice"),
  list(id = "c11", seed = 21, n = 15, N = 300, chg = c("x", "g"), nchg = 80, form = "inertia(net) + recip(net) + indeg(net) + same(nodes$g) + sim(nodes$x) + alter(nodes$x)", sub = "choice"),
  list(id = "c12", seed = 22, n = 15, N = 300, chg = "x", nchg = 60, form = "1 + indeg(net) + ego(nodes$x)", sub = "rate"),
  list(id = "c13", seed = 23, n = 8, N = 80, chg = "x", edge = FALSE, form = "1 + outdeg(net) + ego(nodes$x)", sub = "rate"),
  list(id = "c14", seed = 24, n = 8, N = 80, chg = "x", inc = TRUE, form = "indeg(net) + ego(nodes$x)", sub = "rate"),
  list(id = "c15", seed = 25, n = 7, N = 60, chg = "x", edge = FALSE, form = "inertia(net) + sim(nodes$x) + alter(nodes$x)", sub = "choice")
)
res <- NULL
for (cc in cases) {
  cs <- mkcase(cc$seed, cc$n, cc$N, cc$chg, inc = isTRUE(cc$inc), directed = !isFALSE(cc$directed),
               nchg = if (is.null(cc$nchg)) 15 else cc$nchg, edge = !isFALSE(cc$edge))
  write.csv(cs$ev[, c("time", "sender", "receiver")], file.path(out, paste0("events_", cc$id, ".csv")), row.names = FALSE)
  write.csv(cs$nodes[, c("label", "x", "y", "g")], file.path(out, paste0("nodes_", cc$id, ".csv")), row.names = FALSE)
  for (v in names(cs$chl)) write.csv(cs$chl[[v]], file.path(out, paste0("changes_", v, "_", cc$id, ".csv")), row.names = FALSE)
  f <- tryCatch(fitgf(cs, cc$form, cc$sub), error = function(e) { cat(cc$id, "ERROR", conditionMessage(e), "\n"); NULL })
  if (is.null(f)) next
  cat(cc$id, ":", paste(sprintf("%.8f", c(f$b, f$ll)), collapse = " "), "\n")
  nm <- names(f$b)
  res <- rbind(res, data.frame(case = cc$id, term = c(nm, "ll"), estimate = c(f$b, f$ll), se = c(f$se, NA)))
}
write.csv(res, file.path(out, "results.csv"), row.names = FALSE)
