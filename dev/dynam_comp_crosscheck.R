# nwdynam composition change (actors who join or leave) against goldfish (R 4.x, goldfish 1.6.12).
# Usage: Rscript dev/dynam_comp_crosscheck.R <output directory>
# For every case it writes events_<case>.csv (time, sender, receiver), nodes_<case>.csv (label, present0, x, y, g: the values at the
# start; present0 = presence at the start), comp_<case>.csv (time, node, value: 1 joins, 0 leaves), optionally changes_x_<case>.csv
# (time-varying covariate, replace) and results.csv (case, term, estimate, se, ll) of the goldfish fit; cases.csv holds the nwdynam
# call. The Stata side (cscripts/test_nwdynam_composition.do) reads them.
suppressMessages(library(goldfish))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# presence of node `nd` at time t: all changes at or before t in list order
pres_at <- function(init, comp, lab, t) {
  p <- setNames(init, lab)
  if (nrow(comp)) for (r in which(comp$time <= t)) p[comp$node[r]] <- comp$value[r]
  p
}

mkcase <- function(seed, n, N, pabs = 0.25, ncomp = 12, tvc = FALSE, ties = TRUE, directed = TRUE, edge = NULL) {
  set.seed(seed)
  lab <- paste0("A", 1:n)
  tm <- round(cumsum(runif(N, 0.3, 3)), 2)
  if (ties) for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  tm <- sort(tm)
  init <- runif(n) > pabs
  # composition changes: random times, some exactly at event times, some before the first / after the last event
  ct <- round(runif(ncomp, min(tm) - 2, max(tm) + 2), 2)
  ct[seq(2, ncomp, by = 3)] <- sample(tm, length(seq(2, ncomp, by = 3)))
  ct[1] <- min(tm) - 1; ct[ncomp] <- max(tm) + 1
  ct <- sort(ct)
  comp <- data.frame(time = ct, node = sample(lab, ncomp, replace = TRUE), value = runif(ncomp) < 0.5, stringsAsFactors = FALSE)
  # fixed edge cases: edge(tm, lab) returns the presence at the start and the composition changes
  if (!is.null(edge)) {
    e <- edge(tm, lab)
    init <- e$init
    comp <- data.frame(time = e$time, node = e$node, value = e$value, stringsAsFactors = FALSE)
    comp <- comp[order(comp$time), ]
  }
  # events among present actors only (goldfish stops otherwise); add joins when fewer than 3 actors are present
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) {
    p <- pres_at(init, comp, lab, tm[k])
    while (sum(p) < 3) {
      nd <- sample(lab[!p], 1)
      comp <- rbind(comp, data.frame(time = tm[k], node = nd, value = TRUE, stringsAsFactors = FALSE))
      comp <- comp[order(comp$time), ]
      p <- pres_at(init, comp, lab, tm[k])
    }
    if (k > 4 && k %% 4 == 1 && p[snd[k - 3]] && p[rcv[k - 3]]) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }
    else { s <- sample(lab[p], 1); snd[k] <- s; rcv[k] <- sample(setdiff(lab[p], s), 1) }
  }
  rownames(comp) <- NULL
  ev <- data.frame(time = tm, sender = snd, receiver = rcv, increment = 1, stringsAsFactors = FALSE)
  nodes <- data.frame(label = lab, present = init, x = round(rnorm(n), 2), y = round(runif(n, 0, 2), 2),
                      g = sample(1:3, n, replace = TRUE), stringsAsFactors = FALSE)
  chx <- NULL
  if (tvc) {
    nc <- 12
    ctx <- sort(round(runif(nc, min(tm) - 1, max(tm) + 1), 2)); ctx[2:3] <- sample(tm, 2)
    ctx <- sort(ctx)
    # goldfish refuses a change event of an absent actor as well: draw the actor among those present at that time
    cn <- vapply(ctx, function(t) { p <- pres_at(init, comp, lab, t); if (sum(p) == 0) NA_character_ else sample(lab[p], 1) }, "")
    keep <- !is.na(cn)
    chx <- data.frame(time = ctx[keep], node = cn[keep], value = round(rnorm(sum(keep)), 2), stringsAsFactors = FALSE)
  }
  list(lab = lab, ev = ev, nodes = nodes, comp = comp, chx = chx, directed = directed)
}

fitgf <- function(cs, form, sub) {
  nodes <- defineNodes(cs$nodes)
  assign("ev", cs$ev, globalenv())
  if (nrow(cs$comp) > 0) {
    cc <- data.frame(time = cs$comp$time, node = cs$comp$node, replace = cs$comp$value, stringsAsFactors = FALSE)
    assign("comp", cc, globalenv())
    nodes <- linkEvents(nodes, changeEvent = comp, attribute = "present")
  }
  if (!is.null(cs$chx)) {
    d <- cs$chx; names(d)[3] <- "replace"
    assign("chg_x", d, globalenv())
    nodes <- linkEvents(nodes, changeEvents = chg_x, attribute = "x")
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

# (terms in the order of nwdynam's effect list, so that e(b) and goldfish's coefficients line up)
pieces <- list(
  choice = list(list("inertia(net)", "inertia"), list("recip(net)", "recip"), list("indeg(net)", "indeg"),
    list("same(nodes$%s)", "same(%s)", "gx"), list("diff(nodes$%s)", "diff(%s)", "xy"), list("sim(nodes$%s)", "sim(%s)", "xy"),
    list("alter(nodes$%s)", "alter(%s)", "xy"), list("outdeg(net)", "outdeg"),
    list("trans(net)", "trans"), list("nodeTrans(net)", "nodetrans"), list("egoAlterInt(nodes$x, nodes$y)", "egoalterint(x y)", "", "xy"), list("tertius(net, nodes$%s)", "tertius(%s)", "xy")),
  rate = list(list("indeg(net)", "indeg"), list("outdeg(net)", "outdeg"), list("ego(nodes$%s)", "ego(%s)", "xy"),
    list("tertius(net, nodes$%s)", "tertius(%s)", "xy")),
  choice_coordination = list(list("inertia(net)", "inertia"), list("indeg(net)", "indeg"), list("same(nodes$%s)", "same(%s)", "gx"),
    list("diff(nodes$%s)", "diff(%s)", "xy"), list("sim(nodes$%s)", "sim(%s)", "xy"), list("alter(nodes$%s)", "alter(%s)", "xy"),
    list("tertius(net, nodes$%s)", "tertius(%s)", "xy"), list("egoAlterInt(nodes$x, nodes$y)", "egoalterint(x y)", "", "xy")))

cases <- list(
  list(id = "c01", seed = 1, n = 6, N = 40, form = "inertia(net) + indeg(net)", opts = "inertia indeg", sub = "choice"),
  list(id = "c02", seed = 2, n = 7, N = 50, form = "inertia(net) + same(nodes$g) + alter(nodes$x)", opts = "inertia same(g) alter(x)", sub = "choice"),
  list(id = "c03", seed = 3, n = 7, N = 50, form = "indeg(net) + ego(nodes$x)", opts = "indeg ego(x)", sub = "rate"),
  list(id = "c04", seed = 4, n = 7, N = 50, form = "1 + indeg(net) + ego(nodes$x)", opts = "indeg ego(x) intercept", sub = "rate"),
  list(id = "c05", seed = 5, n = 7, N = 50, form = "inertia(net) + sim(nodes$x) + alter(nodes$x)", opts = "inertia sim(x) alter(x)", sub = "choice_coordination", directed = FALSE),
  list(id = "c06", seed = 6, n = 8, N = 70, tvc = TRUE, form = "1 + indeg(net) + ego(nodes$x)", opts = "indeg ego(x) intercept change(x chg_x)", sub = "rate"),
  list(id = "c07", seed = 7, n = 8, N = 70, tvc = TRUE, form = "inertia(net) + alter(nodes$x)", opts = "inertia alter(x) change(x chg_x)", sub = "choice"),
  list(id = "c08", seed = 8, n = 8, N = 70, pabs = 0.6, form = "1 + indeg(net)", opts = "indeg intercept", sub = "rate"),
  list(id = "c09", seed = 9, n = 9, N = 80, pabs = 0.05, ncomp = 25, form = "inertia(net) + recip(net) + indeg(net)", opts = "inertia recip indeg", sub = "choice"),
  list(id = "c10", seed = 10, n = 8, N = 70, form = "inertia(net, window = 8) + recip(net) + sim(nodes$x)", opts = "inertiawindow(8) recip sim(x)", sub = "choice")
)
# fixed edge cases (tm = the event times, lab = A1..An); the events are drawn among the actors present
edges <- list(
  # A1, A2 leave before the first event and never come back; A3 joins after the last event (never present in the log)
  e01 = function(tm, lab) list(init = c(TRUE, TRUE, FALSE, rep(TRUE, length(lab) - 3)), time = c(min(tm) - 1, min(tm) - 0.5, max(tm) + 1),
        node = c("A1", "A2", "A3"), value = c(FALSE, FALSE, TRUE)),
  # A1 absent during the whole log (no change), A2 present at the start and leaving at the very first event time, joining at the last
  e02 = function(tm, lab) list(init = c(FALSE, rep(TRUE, length(lab) - 1)), time = c(tm[1], tm[length(tm)]), node = c("A2", "A2"), value = c(FALSE, TRUE)),
  # joins and leaves exactly at event times (the join takes effect for that very event, the leave too), several changes of one actor at one time
  e03 = function(tm, lab) list(init = c(FALSE, TRUE, TRUE, TRUE, FALSE, rep(TRUE, length(lab) - 5)),
        time = c(tm[5], tm[9], tm[9], tm[14], tm[14], tm[20], tm[25]), node = c("A1", "A2", "A2", "A5", "A5", "A1", "A3"),
        value = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)),
  # many actors join late: only 3 present at the start
  e04 = function(tm, lab) list(init = c(rep(TRUE, 3), rep(FALSE, length(lab) - 3)),
        time = tm[c(4, 8, 12, 16, 20, 24)], node = lab[4:9], value = rep(TRUE, 6)),
  # everybody leaves once and returns; leaves and joins in the middle of tied event times
  e05 = function(tm, lab) list(init = rep(TRUE, length(lab)), time = c(tm[7] - 0.01, tm[7], tm[7] + 0.01, tm[30], tm[31]), node = c("A4", "A5", "A4", "A6", "A6"),
        value = c(FALSE, FALSE, TRUE, FALSE, TRUE))
)
fixed <- list(
  list(id = "e01c", seed = 31, n = 8, N = 45, edge = "e01", form = "inertia(net) + indeg(net)", opts = "inertia indeg", sub = "choice"),
  list(id = "e01r", seed = 31, n = 8, N = 45, edge = "e01", form = "indeg(net) + ego(nodes$x)", opts = "indeg ego(x)", sub = "rate"),
  list(id = "e01i", seed = 31, n = 8, N = 45, edge = "e01", form = "1 + indeg(net)", opts = "indeg intercept", sub = "rate"),
  list(id = "e02c", seed = 32, n = 7, N = 45, edge = "e02", form = "inertia(net) + recip(net)", opts = "inertia recip", sub = "choice"),
  list(id = "e02i", seed = 32, n = 7, N = 45, edge = "e02", form = "1 + outdeg(net)", opts = "outdeg intercept", sub = "rate"),
  list(id = "e03c", seed = 33, n = 8, N = 50, edge = "e03", form = "inertia(net) + indeg(net) + alter(nodes$x)", opts = "inertia indeg alter(x)", sub = "choice"),
  list(id = "e03i", seed = 33, n = 8, N = 50, edge = "e03", form = "1 + indeg(net) + ego(nodes$x)", opts = "indeg ego(x) intercept", sub = "rate"),
  list(id = "e03r", seed = 33, n = 8, N = 50, edge = "e03", form = "indeg(net) + ego(nodes$x)", opts = "indeg ego(x)", sub = "rate"),
  list(id = "e03m", seed = 33, n = 8, N = 50, edge = "e03", form = "inertia(net) + sim(nodes$x)", opts = "inertia sim(x)", sub = "choice_coordination", directed = FALSE),
  list(id = "e04c", seed = 34, n = 9, N = 50, edge = "e04", form = "inertia(net) + indeg(net)", opts = "inertia indeg", sub = "choice"),
  list(id = "e04i", seed = 34, n = 9, N = 50, edge = "e04", form = "1 + indeg(net)", opts = "indeg intercept", sub = "rate"),
  list(id = "e05c", seed = 35, n = 8, N = 50, edge = "e05", form = "inertia(net) + recip(net)", opts = "inertia recip", sub = "choice"),
  list(id = "e05i", seed = 35, n = 8, N = 50, edge = "e05", form = "1 + indeg(net)", opts = "indeg intercept", sub = "rate"),
  list(id = "e05m", seed = 35, n = 8, N = 50, edge = "e05", form = "inertia(net) + alter(nodes$x)", opts = "inertia alter(x)", sub = "choice_coordination", directed = FALSE)
)
cases <- c(cases, fixed)
set.seed(2026)
nrand <- 50
for (k in 1:nrand) {
  sub <- c("choice", "rate", "choice", "rate", "choice_coordination", "choice")[(k - 1) %% 6 + 1]
  icpt <- (sub == "rate" && k %% 4 == 2)
  pc <- pieces[[sub]]
  take <- which(runif(length(pc)) < 0.4)
  if (length(take) == 0) take <- sample(length(pc), 1)
  take <- take[seq_len(min(length(take), 4))]
  forms <- character(); opts <- character(); used <- character()
  for (t in take) {
    pp <- pc[[t]]
    if (length(pp) == 2) { forms <- c(forms, pp[[1]]); opts <- c(opts, pp[[2]]) }
    else if (length(pp) == 3) {
      v <- sample(strsplit(pp[[3]], "")[[1]], 1)
      forms <- c(forms, sprintf(pp[[1]], v)); opts <- c(opts, sprintf(pp[[2]], v)); used <- c(used, v)
    } else { forms <- c(forms, pp[[1]]); opts <- c(opts, pp[[2]]); used <- c(used, "x", "y") }
  }
  tvc <- ("x" %in% used) && (k %% 5 == 0)
  if (tvc) opts <- c(opts, "change(x chg_x)")
  cases[[length(cases) + 1]] <- list(id = sprintf("r%02d", k), seed = 700 + k, n = sample(6:12, 1), N = sample(40:110, 1), tvc = tvc,
    pabs = sample(c(0.1, 0.25, 0.45), 1), ncomp = sample(6:30, 1),
    form = paste(c(if (icpt) "1", forms), collapse = " + "), sub = sub, directed = (sub != "choice_coordination"),
    opts = paste(c(opts, if (icpt) "intercept"), collapse = " "))
}

res <- NULL; cdf <- NULL
for (cc in cases) {
  cs <- mkcase(cc$seed, cc$n, cc$N, pabs = if (is.null(cc$pabs)) 0.25 else cc$pabs, ncomp = if (is.null(cc$ncomp)) 12 else cc$ncomp,
               tvc = isTRUE(cc$tvc), directed = !isFALSE(cc$directed), edge = if (is.null(cc$edge)) NULL else edges[[cc$edge]])
  f <- tryCatch(fitgf(cs, cc$form, cc$sub), error = function(e) { cat(cc$id, "SKIPPED", conditionMessage(e), "\n"); NULL })
  if (is.null(f) || any(!is.finite(f$se)) || max(abs(f$b[names(f$b) != "Intercept"])) > 8 || max(f$se) > 3) { if (!is.null(f)) cat(cc$id, "SKIPPED (unstable)\n"); next }
  write.csv(cs$ev[, c("time", "sender", "receiver")], file.path(out, paste0("events_", cc$id, ".csv")), row.names = FALSE)
  nd <- cs$nodes[, c("label", "x", "y", "g")]; nd$present0 <- as.integer(cs$nodes$present)
  write.csv(nd, file.path(out, paste0("nodes_", cc$id, ".csv")), row.names = FALSE)
  write.csv(data.frame(time = cs$comp$time, node = cs$comp$node, value = as.integer(cs$comp$value)), file.path(out, paste0("comp_", cc$id, ".csv")), row.names = FALSE)
  if (!is.null(cs$chx)) write.csv(cs$chx, file.path(out, paste0("changes_x_", cc$id, ".csv")), row.names = FALSE)
  cat(cc$id, ":", paste(sprintf("%.8f", c(f$b, f$ll)), collapse = " "), "\n")
  nm <- names(f$b)
  res <- rbind(res, data.frame(case = cc$id, term = c(nm, "ll"), estimate = c(f$b, f$ll), se = c(f$se, NA)))
  cdf <- rbind(cdf, data.frame(id = cc$id, sub = cc$sub, opts = cc$opts, directed = as.integer(!isFALSE(cc$directed)),
    tolb = if (cc$sub == "choice_coordination") 1e-5 else 1e-6, tolse = if (cc$sub == "choice_coordination") 0.1 else 1e-6))
}
write.csv(res, file.path(out, "results.csv"), row.names = FALSE)
write.csv(cdf, file.path(out, "cases.csv"), row.names = FALSE)
