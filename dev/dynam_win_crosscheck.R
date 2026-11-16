# nwdynam windows, weights, two-mode and opportunities combinations against goldfish (R 4.x, goldfish 1.6.12).
# Usage: Rscript dev/dynam_win_crosscheck.R <output directory> [sections, e.g. a b c d e f]
# Sections: a = windowed rate effects with intercept, b = windows on tertius, c = choice_coordination weighted/windowed,
# d = two-mode attribute and closure effects, e = opportunities with attribute effects, f = time-varying covariates in two-mode
# and with opportunities. For every case it writes events_<case>.csv (time, sender, receiver, value), nodes_<case>.csv (label, x, y, g),
# changes_<var>_<case>.csv (time, node, value), opp_<case>.csv (event, actor) and, per section, results_<sec>.csv (case, term,
# estimate, se; the last row of a case is "ll") and cases_<sec>.csv (case, sub, options of nwdynam, kind, tolerances).
# The Stata side (dev/dynam_win_crosscheck.do) reads them.
suppressMessages(library(goldfish))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "."
secs <- if (length(args) >= 2) args[-1] else c("a", "b", "c", "d", "e", "f")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- one random log. kind: "one" (directed), "und" (undirected), "two" (two-mode: n1 senders, n2 receivers)
mkcase <- function(seed, n = 8, N = 60, kind = "one", n1 = 6, n2 = 5, weights = FALSE, inttime = FALSE, ties = TRUE,
                   chg = character(0), nchg = 15, inc = FALSE, opp = FALSE, oppsize = 0.6, canon = FALSE) {
  set.seed(seed)
  if (kind == "two") {
    l1 <- paste0("P", 1:n1); l2 <- paste0("O", 1:n2); lab <- c(l1, l2); n <- n1 + n2
  } else { lab <- paste0("A", 1:n); l1 <- lab; l2 <- lab }
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) {
    if (kind == "two") { snd[k] <- sample(l1, 1); rcv[k] <- sample(l2, 1) }
    else { s <- sample(lab, 1); snd[k] <- s; rcv[k] <- sample(setdiff(lab, s), 1) }
  }
  for (k in seq(5, N, by = 4)) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }   # repetitions
  if (inttime) { tm <- cumsum(sample(0:3, N, TRUE, prob = c(.15, .35, .3, .2))); tm <- tm + 1 }
  else {
    tm <- round(cumsum(runif(N, 0.3, 3)) * 4) / 4   # multiples of 0.25: the event times of the core are float, so only such times are exact
    if (ties) for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  }
  tm <- sort(tm)
  val <- if (weights) sample(1:3, N, TRUE, prob = c(.6, .25, .15)) else rep(1, N)
  if (kind == "und" && canon) { # undirected, every pair in one orientation only (lower index sends)
    sw <- match(snd, lab) > match(rcv, lab); tmp <- snd[sw]; snd[sw] <- rcv[sw]; rcv[sw] <- tmp
  }
  if (kind == "und") { # undirected: orientation arbitrary unless canon
    ev <- data.frame(time = tm, sender = snd, receiver = rcv, increment = val, stringsAsFactors = FALSE)
  } else ev <- data.frame(time = tm, sender = snd, receiver = rcv, increment = val, stringsAsFactors = FALSE)
  nodes <- data.frame(label = lab, present = TRUE, x = round(rnorm(n), 2), y = round(runif(n, 0, 2), 2),
                      g = sample(1:3, n, replace = TRUE), stringsAsFactors = FALSE)
  chl <- list()
  for (v in chg) {
    ct <- sort(round(runif(nchg, min(tm) - 2, max(tm) + 2) * 4) / 4)
    ct[2:4] <- tm[c(3, 10, 20)]; ct[5:6] <- tm[15]; ct <- sort(ct)
    pool <- lab
    cn <- sample(pool, nchg, replace = TRUE)
    if (v == "g") cv <- sample(1:3, nchg, replace = TRUE) else cv <- round(rnorm(nchg), 2)
    if (inc && v != "g") cv <- round(rnorm(nchg, 0, 0.5), 2)
    chl[[v]] <- data.frame(time = ct, node = cn, value = cv, stringsAsFactors = FALSE)
  }
  oppdf <- NULL
  if (opp) {
    rows <- list()
    for (k in 1:N) {
      cand <- if (kind == "two") match(l2, lab) else setdiff(seq_len(n), match(snd[k], lab))
      rr <- match(rcv[k], lab)
      av <- unique(c(rr, cand[runif(length(cand)) < oppsize]))
      rows[[k]] <- data.frame(event = k, actor = av)
    }
    oppdf <- do.call(rbind, rows)
  }
  list(lab = lab, l1 = l1, l2 = l2, n1 = n1, n2 = n2, kind = kind, ev = ev, nodes = nodes, chl = chl, inc = inc, opp = oppdf)
}

# ---- fit in goldfish. form: formula with objects net, nd (one-mode nodes) or nd1, nd2 (two-mode), dep is the dependent object
fitgf <- function(cs, form, sub, intercept = FALSE) {
  ev <- cs$ev
  assign("ev", ev, globalenv())
  mkn <- function(df, labs, extra = NULL) {
    d <- df[df$label %in% labs, ]
    nd <- defineNodes(d)
    for (v in names(cs$chl)) {
      ch <- cs$chl[[v]]; ch <- ch[ch$node %in% labs, ]
      if (nrow(ch) == 0) next
      names(ch)[3] <- if (cs$inc) "increment" else "replace"
      nm <- paste0("chg_", v, "_", labs[1])
      assign(nm, ch, globalenv())
      nd <- eval(substitute(linkEvents(nd, changeEvents = CH, attribute = V), list(CH = as.name(nm), V = v)))
    }
    nd
  }
  if (cs$kind == "two") {
    nd1 <- mkn(cs$nodes, cs$l1); nd2 <- mkn(cs$nodes, cs$l2)
    assign("nd1", nd1, globalenv()); assign("nd2", nd2, globalenv())
    net <- defineNetwork(nodes = nd1, nodes2 = nd2, directed = TRUE)
    net <- linkEvents(net, changeEvent = ev, nodes = nd1, nodes2 = nd2)
    dep <- defineDependentEvents(events = ev, nodes = nd1, nodes2 = nd2, defaultNetwork = net)
  } else {
    nd <- mkn(cs$nodes, cs$lab)
    assign("nd", nd, globalenv())
    net <- defineNetwork(nodes = nd, directed = (cs$kind == "one"))
    net <- linkEvents(net, changeEvent = ev, nodes = nd)
    dep <- defineDependentEvents(events = ev, nodes = nd, defaultNetwork = net)
  }
  assign("net", net, globalenv()); assign("dep", dep, globalenv())
  ei <- list(maxIterations = 200, initialDamping = 10, maxScoreStopCriterion = 1e-9)
  if (!is.null(cs$opp)) {
    # opportunities: a list with, for every event (in order), the available actors (indices)
    ol <- split(cs$opp$actor, cs$opp$event)
    ol <- lapply(1:nrow(ev), function(k) as.numeric(ol[[as.character(k)]]))
    ei$opportunitiesList <- ol
  }
  m <- suppressWarnings(estimate(as.formula(paste("dep ~", form), env = globalenv()), model = "DyNAM", subModel = sub, estimationInit = ei))
  list(b = coef(m), se = m$standardErrors, ll = as.numeric(m$logLikelihood))
}

res <- list(); cdf <- list()
# run a list of case specs; each spec: sec, id, sub, form, opts, plus mkcase arguments in `mk`
runspec <- function(sp) {
  cs <- do.call(mkcase, c(list(seed = sp$seed), sp$mk))
  f <- tryCatch(fitgf(cs, sp$form, sp$sub), error = function(e) { cat(sp$id, "SKIPPED:", conditionMessage(e), "\n"); NULL })
  if (is.null(f)) return(invisible())
  if (any(!is.finite(f$se)) || max(abs(f$b[names(f$b) != "Intercept"])) > 8 || max(f$se) > (if (sp$sub == "choice_coordination") 0.6 else 3)) { cat(sp$id, "SKIPPED (unstable)\n"); return(invisible()) }
  id <- sp$id
  write.csv(cs$ev[, c("time", "sender", "receiver", "increment")], file.path(out, paste0("events_", id, ".csv")), row.names = FALSE)
  write.csv(cs$nodes[, c("label", "x", "y", "g")], file.path(out, paste0("nodes_", id, ".csv")), row.names = FALSE)
  for (v in names(cs$chl)) write.csv(cs$chl[[v]], file.path(out, paste0("changes_", v, "_", id, ".csv")), row.names = FALSE)
  if (!is.null(cs$opp)) write.csv(cs$opp, file.path(out, paste0("opp_", id, ".csv")), row.names = FALSE)
  cat(id, ":", paste(sprintf("%.8f", c(f$b, f$ll)), collapse = " "), "\n")
  coord <- (sp$sub == "choice_coordination")
  res[[length(res) + 1]] <<- data.frame(sec = sp$sec, case = id, term = c(names(f$b), "ll"), estimate = c(f$b, f$ll), se = c(f$se, NA))
  cdf[[length(cdf) + 1]] <<- data.frame(sec = sp$sec, id = id, sub = sp$sub, opts = sp$opts, kind = cs$kind,
    chg = paste(names(cs$chl), collapse = " "), opp = as.integer(!is.null(cs$opp)), inc = as.integer(isTRUE(cs$inc)),
    tolb = if (coord) 5e-4 else 1e-6, tolse = if (coord) 0.1 else 1e-6)
}

# ---- random specs from pieces: a piece = list(goldfish term, nwdynam option[, variable letters])
randspecs <- function(sec, pieces, nrand, seed0, sub, mkargs, kind = "one", maxeff = 3, icpt = FALSE, idoff = 0, groups = NULL) {
  set.seed(seed0)
  sp <- list()
  for (k in 1:nrand) {
    take <- which(runif(length(pieces)) < 0.5)
    if (length(take) == 0) take <- sample(length(pieces), 1)
    if (!is.null(groups)) take <- take[!duplicated(groups[take])]
    take <- take[seq_len(min(length(take), maxeff))]
    forms <- character(); opts <- character()
    for (t in take) {
      pp <- pieces[[t]]
      if (length(pp) >= 3) {
        v <- sample(strsplit(pp[[3]], "")[[1]], 1)
        forms <- c(forms, gsub("%s", v, pp[[1]], fixed = TRUE)); opts <- c(opts, gsub("%s", v, pp[[2]], fixed = TRUE))
      } else { forms <- c(forms, pp[[1]]); opts <- c(opts, pp[[2]]) }
    }
    ic <- icpt && (k %% 2 == 0)
    mk <- mkargs(k)
    fm <- paste(c(if (ic) "1", forms), collapse = " + "); op <- paste(c(opts, if (ic) "intercept"), collapse = " ")
    while (grepl("%W", fm, fixed = TRUE)) { w <- sample(c(2, 3, 4.5, 6, 9), 1); fm <- sub("%W", w, fm, fixed = TRUE); op <- sub("%W", w, op, fixed = TRUE) }
    if (length(mk$chg) > 0) mk$chg <- Filter(function(v) grepl(paste0("[( ]", v, "[) ]|[( ]", v, "\\)|\\(", v, "$"), op), mk$chg)
    sp[[k]] <- list(sec = sec, id = sprintf("%s%02d", sec, k + idoff), seed = seed0 * 100 + k, sub = sub, mk = mk, form = fm, opts = op)
  }
  sp
}

specs <- list()
# ================= section a: windowed rate effects together with the intercept =================
if ("a" %in% secs) {
  pa <- list(list("indeg(net, window = %W)", "indegwindow(%W)"), list("outdeg(net, window = %W)", "outdegwindow(%W)"),
             list("ego(nd$%s)", "ego(%s)", "xy"))
  sa <- list(
    list(sec = "a", id = "a01", seed = 1, sub = "rate", mk = list(n = 7, N = 50), form = "1 + indeg(net, window = 6)", opts = "indegwindow(6) intercept"),
    list(sec = "a", id = "a02", seed = 2, sub = "rate", mk = list(n = 7, N = 50), form = "1 + outdeg(net, window = 6)", opts = "outdegwindow(6) intercept"),
    list(sec = "a", id = "a03", seed = 3, sub = "rate", mk = list(n = 8, N = 70, inttime = TRUE), form = "1 + indeg(net, window = 5) + outdeg(net, window = 4)", opts = "indegwindow(5) outdegwindow(4) intercept"),
    list(sec = "a", id = "a04", seed = 4, sub = "rate", mk = list(n = 8, N = 70, weights = TRUE), form = "1 + indeg(net, window = 6, weighted = TRUE) + ego(nd$x)", opts = "indegwindow(6) weightedindeg ego(x) intercept"),
    list(sec = "a", id = "a05", seed = 5, sub = "rate", mk = list(n = 8, N = 70, weights = TRUE, inttime = TRUE), form = "1 + outdeg(net, window = 3, weighted = TRUE)", opts = "outdegwindow(3) weightedoutdeg intercept"))
  # random: windows on indeg/outdeg with ego, weights, integer times, ties, change events (also after the last event)
  sa2 <- randspecs("a", pa, 30, 11, "rate", function(k) list(n = sample(6:10, 1), N = sample(40:90, 1), weights = (k %% 3 == 0), inttime = (k %% 4 == 0),
         chg = if (k %% 3 == 1) "x" else character(0), nchg = 12, inc = (k %% 6 == 4)), icpt = TRUE, idoff = 10)
  # all with the intercept
  for (i in seq_along(sa2)) if (!grepl("intercept", sa2[[i]]$opts)) { sa2[[i]]$form <- paste("1 +", sa2[[i]]$form); sa2[[i]]$opts <- paste(sa2[[i]]$opts, "intercept") }
  # a window shorter than the gap to the changes after the last event, and changes with ties between events
  sa3 <- list(
    list(sec = "a", id = "a60", seed = 60, sub = "rate", mk = list(n = 7, N = 50, chg = "x", nchg = 14), form = "1 + indeg(net, window = 0.75) + ego(nd$x)", opts = "indegwindow(0.75) ego(x) intercept"),
    list(sec = "a", id = "a61", seed = 61, sub = "rate", mk = list(n = 7, N = 50, chg = "x", nchg = 14, inttime = TRUE), form = "1 + outdeg(net, window = 1) + ego(nd$x)", opts = "outdegwindow(1) ego(x) intercept"),
    list(sec = "a", id = "a62", seed = 62, sub = "rate", mk = list(n = 9, N = 80, chg = "y", nchg = 20, weights = TRUE), form = "1 + indeg(net, window = 1.25, weighted = TRUE) + outdeg(net, window = 2.5) + ego(nd$y)", opts = "indegwindow(1.25) weightedindeg outdegwindow(2.5) ego(y) intercept"))
  specs <- c(specs, sa, sa2, sa3)
}
# ================= section c: choice_coordination, windowed and weighted inertia / indeg =================
if ("c" %in% secs) {
  pc <- list(list("inertia(net, window = %W)", "inertiawindow(%W)"), list("inertia(net, weighted = TRUE)", "weightedinertia"),
             list("indeg(net, window = %W)", "indegwindow(%W)"), list("indeg(net, weighted = TRUE)", "weightedindeg"),
             list("inertia(net, window = %W, weighted = TRUE)", "inertiawindow(%W) weightedinertia"),
             list("indeg(net, window = %W, weighted = TRUE)", "indegwindow(%W) weightedindeg"),
             list("same(nd$g)", "same(g)"), list("sim(nd$%s)", "sim(%s)", "xy"), list("alter(nd$%s)", "alter(%s)", "xy"),
             list("trans(net)", "trans"), list("tertius(net, nd$%s)", "tertius(%s)", "xy"))
  sc <- list(
    list(sec = "c", id = "c01", seed = 1, sub = "choice_coordination", mk = list(kind = "und", n = 7, N = 50), form = "inertia(net, window = 8)", opts = "inertiawindow(8)"),
    list(sec = "c", id = "c02", seed = 2, sub = "choice_coordination", mk = list(kind = "und", n = 7, N = 50, weights = TRUE), form = "inertia(net, weighted = TRUE)", opts = "weightedinertia"),
    list(sec = "c", id = "c03", seed = 3, sub = "choice_coordination", mk = list(kind = "und", n = 7, N = 50), form = "indeg(net, window = 8)", opts = "indegwindow(8)"),
    list(sec = "c", id = "c04", seed = 4, sub = "choice_coordination", mk = list(kind = "und", n = 7, N = 50, weights = TRUE), form = "indeg(net, weighted = TRUE)", opts = "weightedindeg"))
  sc2 <- randspecs("c", pc, 40, 31, "choice_coordination", function(k) list(kind = "und", n = sample(6:10, 1), N = sample(40:90, 1),
         weights = (k %% 2 == 0), inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 10, groups = c(1, 1, 2, 2, 1, 2, 3, 4, 5, 6, 7))
  sc3 <- list(
    list(sec = "c", id = "c61", seed = 61, sub = "choice_coordination", mk = list(kind = "und", n = 8, N = 70), form = "inertia(net, window = 6) + indeg(net, window = 4.5)", opts = "inertiawindow(6) indegwindow(4.5)"),
    list(sec = "c", id = "c62", seed = 62, sub = "choice_coordination", mk = list(kind = "und", n = 8, N = 70), form = "inertia(net, window = 6) + indeg(net)", opts = "inertiawindow(6) indeg"),
    list(sec = "c", id = "c63", seed = 63, sub = "choice_coordination", mk = list(kind = "und", n = 8, N = 70), form = "inertia(net) + indeg(net, window = 4.5)", opts = "inertia indegwindow(4.5)"),
    list(sec = "c", id = "c64", seed = 64, sub = "choice_coordination", mk = list(kind = "und", n = 8, N = 70), form = "inertia(net, window = 6) + sim(nd$x)", opts = "inertiawindow(6) sim(x)"),
    list(sec = "c", id = "c65", seed = 65, sub = "choice_coordination", mk = list(kind = "und", n = 8, N = 70), form = "indeg(net, window = 4.5) + sim(nd$x)", opts = "indegwindow(4.5) sim(x)"))
  sc4 <- randspecs("c", pc, 40, 41, "choice_coordination", function(k) list(kind = "und", n = sample(6:10, 1), N = sample(40:90, 1), canon = TRUE,
         weights = (k %% 2 == 0), inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 100, groups = c(1, 1, 2, 2, 1, 2, 3, 4, 5, 6, 7))
  specs <- c(specs, sc, sc2, sc3, sc4)
}
for (sp in specs) runspec(sp)
for (sec in secs) {
  r <- do.call(rbind, Filter(function(d) d$sec[1] == sec, res)); cd <- do.call(rbind, Filter(function(d) d$sec[1] == sec, cdf))
  if (!is.null(r)) { write.csv(r[, -1], file.path(out, paste0("results_", sec, ".csv")), row.names = FALSE); write.csv(cd[, -1], file.path(out, paste0("cases_", sec, ".csv")), row.names = FALSE) }
}
