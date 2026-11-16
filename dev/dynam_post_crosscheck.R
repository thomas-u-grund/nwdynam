# nwdynam postestimation against goldfish (R 4.x, goldfish 1.6.12, changepoint 2.3): the log likelihood of every event
# (estimation option returnIntervalLogL), the event probabilities (returnEventProbabilities), logLik / AIC / BIC, examineOutliers and
# examineChangepoints (the package changepoint as examineChangepoints calls it).
# Usage: Rscript dev/dynam_post_crosscheck.R <output directory>      (the tests read the files; they do not need R)
# Families of random event logs (every case is fitted in goldfish and written to files):
#   r*  one-mode logs for the choice, rate, rate-with-intercept and choice_coordination sub-models with composition change (actors who join
#       or leave), time-varying covariates and windows in random combinations (also two-mode: m*)
#   o*  choice with an opportunity list per event (oppframe())
#   x*  netterms(): effects on other networks (tie, trans, cycle, four, mixed closures, windows)
# Files per case <id>: events_<id>.csv (time, sender, receiver, increment), nodes_<id>.csv (label, x, y, g, present0), comp_<id>.csv,
# changes_x_<id>.csv, avail_<id>.csv (opportunities), exo_<id>_<name>.csv; reference values: ref_ll_<id>.csv (one row per record of
# goldfish's intervalLogL: type 1 = dependent event, 2 = right-censored interval), ref_pr_<id>.csv (event, actor / actor1, actor2, value),
# and the files post_ic.csv, post_outliers.csv, post_cpts.csv, results.csv, cases.csv.
suppressMessages(library(goldfish)); suppressMessages(library(changepoint))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
EI <- list(maxIterations = 200, initialDamping = 10, maxScoreStopCriterion = 1e-9)

pres_at <- function(init, comp, lab, t) {
  p <- setNames(init, lab)
  if (nrow(comp)) for (r in which(comp$time <= t)) p[comp$node[r]] <- comp$value[r]
  p
}

# ---------------------------------------------------------------- random one-mode logs (composition change, time-varying covariate)
mkcase <- function(seed, n, N, pabs = 0.25, ncomp = 12, tvc = FALSE, ties = TRUE, directed = TRUE, nocomp = FALSE, cycle = FALSE) {
  set.seed(seed)
  lab <- paste0("A", 1:n)
  tm <- round(cumsum(runif(N, 0.3, 3)), 2)
  if (ties) for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  tm <- sort(tm)
  init <- if (nocomp) rep(TRUE, n) else runif(n) > pabs
  if (nocomp) comp <- data.frame(time = numeric(0), node = character(0), value = logical(0), stringsAsFactors = FALSE)
  else {
    ct <- round(runif(ncomp, min(tm) - 2, max(tm) + 2), 2)
    ct[seq(2, ncomp, by = 3)] <- sample(tm, length(seq(2, ncomp, by = 3)))
    ct[1] <- min(tm) - 1; ct[ncomp] <- max(tm) + 1
    ct <- sort(ct)
    comp <- data.frame(time = ct, node = sample(lab, ncomp, replace = TRUE), value = runif(ncomp) < 0.5, stringsAsFactors = FALSE)
  }
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) {
    p <- pres_at(init, comp, lab, tm[k])
    while (sum(p) < 3) {
      nd <- sample(lab[!p], 1)
      comp <- rbind(comp, data.frame(time = tm[k], node = nd, value = TRUE, stringsAsFactors = FALSE))
      comp <- comp[order(comp$time), ]
      p <- pres_at(init, comp, lab, tm[k])
    }
    if (cycle && k <= n) { snd[k] <- lab[k]; rcv[k] <- lab[if (k == n) 1 else k + 1] }
    else if (k > 4 && k %% 4 == 1 && p[snd[k - 3]] && p[rcv[k - 3]]) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }
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
    cn <- vapply(ctx, function(t) { p <- pres_at(init, comp, lab, t); if (sum(p) == 0) NA_character_ else sample(lab[p], 1) }, "")
    keep <- !is.na(cn)
    chx <- data.frame(time = ctx[keep], node = cn[keep], value = round(rnorm(sum(keep)), 2), stringsAsFactors = FALSE)
  }
  list(lab = lab, ev = ev, nodes = nodes, comp = comp, chx = chx, directed = directed, n = n, N = N, nocomp = nocomp)
}

# goldfish fit with the post-estimation results; returns the model and the pieces needed for the references
fitgf <- function(cs, form, sub, opp = NULL) {
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
  for (nm in names(cs$exo)) {
    e <- cs$exo[[nm]]
    if (e$type == "stat") x <- defineNetwork(matrix = e$M, nodes = nodes, directed = e$directed)
    else {
      evn <- paste0("exoev_", nm); assign(evn, e$ev, globalenv())
      x <- defineNetwork(nodes = nodes, directed = e$directed)
      x <- eval(substitute(linkEvents(x, changeEvent = EV, nodes = nodes), list(EV = as.name(evn))))
    }
    assign(nm, x, globalenv())
  }
  dep <- defineDependentEvents(events = ev, nodes = nodes, defaultNetwork = net)
  assign("dep", dep, globalenv())
  fm <- as.formula(paste("dep ~", form), env = globalenv())
  ei <- c(EI, list(returnIntervalLogL = TRUE, returnEventProbabilities = TRUE))
  if (!is.null(opp)) ei$opportunitiesList <- opp
  m <- suppressWarnings(do.call(estimate, list(fm, model = "DyNAM", subModel = sub, estimationInit = ei)))
  pi <- list(maxIterations = 0)
  prep <- suppressWarnings(do.call(estimate, list(fm, model = "DyNAM", subModel = sub, preprocessingOnly = TRUE)))
  list(m = m, prep = prep, b = coef(m), se = m$standardErrors, ll = as.numeric(m$logLikelihood), dep = dep, nodes = nodes)
}

fit2 <- function(cs, form, sub) {
  n1 <- length(cs$l1)
  nd1 <- defineNodes(data.frame(label = cs$l1, present = cs$init[1:n1], stringsAsFactors = FALSE))
  nd2 <- defineNodes(data.frame(label = cs$l2, present = cs$init[-(1:n1)], stringsAsFactors = FALSE))
  c1 <- cs$comp[cs$comp$node %in% cs$l1, ]; c2 <- cs$comp[cs$comp$node %in% cs$l2, ]
  if (nrow(c1)) { assign("comp1", data.frame(time = c1$time, node = c1$node, replace = c1$value, stringsAsFactors = FALSE), globalenv()); nd1 <- linkEvents(nd1, changeEvent = comp1, attribute = "present") }
  if (nrow(c2)) { assign("comp2", data.frame(time = c2$time, node = c2$node, replace = c2$value, stringsAsFactors = FALSE), globalenv()); nd2 <- linkEvents(nd2, changeEvent = comp2, attribute = "present") }
  assign("ev", cs$ev, globalenv()); assign("nd1", nd1, globalenv()); assign("nd2", nd2, globalenv())
  net <- defineNetwork(nodes = nd1, nodes2 = nd2, directed = TRUE)
  net <- linkEvents(net, changeEvent = ev, nodes = nd1, nodes2 = nd2); assign("net", net, globalenv())
  dep <- defineDependentEvents(events = ev, nodes = nd1, nodes2 = nd2, defaultNetwork = net); assign("dep", dep, globalenv())
  fm <- as.formula(paste("dep ~", form), env = globalenv())
  ei <- c(EI, list(returnIntervalLogL = TRUE, returnEventProbabilities = TRUE))
  m <- suppressWarnings(do.call(estimate, list(fm, model = "DyNAM", subModel = sub, estimationInit = ei)))
  prep <- suppressWarnings(do.call(estimate, list(fm, model = "DyNAM", subModel = sub, preprocessingOnly = TRUE)))
  list(m = m, prep = prep, b = coef(m), se = m$standardErrors, ll = as.numeric(m$logLikelihood), dep = dep)
}

mk2 <- function(seed, n1, n2, N, pabs = 0.25, ncomp = 10) {
  set.seed(seed)
  l1 <- paste0("P", 1:n1); l2 <- paste0("O", 1:n2); lab <- c(l1, l2)
  tm <- round(cumsum(runif(N, 0.3, 3)), 2)
  for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  init <- runif(n1 + n2) > pabs
  ct <- sort(round(runif(ncomp, min(tm) - 2, max(tm) + 2), 2))
  ct[seq(2, ncomp, by = 3)] <- sort(sample(tm, length(seq(2, ncomp, by = 3))))
  ct <- sort(ct)
  comp <- data.frame(time = ct, node = sample(lab, ncomp, replace = TRUE), value = runif(ncomp) < 0.5, stringsAsFactors = FALSE)
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) {
    p <- pres_at(init, comp, lab, tm[k])
    while (sum(p[l1]) < 1 || sum(p[l2]) < 2) {
      nd <- sample(lab[!p], 1)
      comp <- rbind(comp, data.frame(time = tm[k], node = nd, value = TRUE, stringsAsFactors = FALSE)); comp <- comp[order(comp$time), ]
      p <- pres_at(init, comp, lab, tm[k])
    }
    snd[k] <- sample(l1[p[l1]], 1); rcv[k] <- sample(l2[p[l2]], 1)
  }
  used <- unique(c(snd, rcv)); keep <- lab %in% used
  l1 <- l1[l1 %in% used]; l2 <- l2[l2 %in% used]; init <- init[keep]; lab <- lab[keep]
  comp <- comp[comp$node %in% used, ]
  rownames(comp) <- NULL
  list(l1 = l1, l2 = l2, lab = lab, init = init, comp = comp,
       ev = data.frame(time = tm, sender = snd, receiver = rcv, increment = 1, stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------- references from a goldfish fit
sortpair <- function(a, b) { o <- order(c(a, b), method = "radix"); c(a, b)[o] }

# the candidate labels of event k as goldfish returns them (vector of the probabilities) -> data frame (event, actor, value)
refprobs <- function(f, cs, sub, twomode = FALSE, oppmat = NULL) {
  m <- f$m; prep <- f$prep
  dep <- which(prep$orderEvents == 1)
  P <- m$eventProbabilities
  tm <- prep$eventTime[dep]
  snd <- prep$eventSender[dep]; rcv <- prep$eventReceiver[dep]
  rows <- list()
  if (twomode) { lab1 <- cs$l1; lab2 <- cs$l2; init <- cs$init; labs <- cs$lab; comp <- cs$comp }
  else { lab1 <- cs$lab; lab2 <- cs$lab; init <- cs$nodes$present; labs <- cs$lab; comp <- cs$comp }
  n <- length(lab1)
  for (k in seq_along(dep)) {
    v <- P[[dep[k]]]
    t <- prep$eventTime[dep[k]]
    p <- pres_at(c(init), comp, labs, t)
    s <- snd[k]
    if (sub == "choice_coordination") {
      keep <- if (nrow(comp) > 0 && !isTRUE(cs$nocomp)) p[lab1] else rep(TRUE, n)
      L <- lab1[keep]
      M <- v
      for (a in seq_len(length(L) - 1)) for (b in (a + 1):length(L))
        rows[[length(rows) + 1]] <- data.frame(event = k, actor1 = sortpair(L[a], L[b])[1], actor2 = sortpair(L[a], L[b])[2], value = M[a, b], stringsAsFactors = FALSE)
    } else if (sub == "choice") {
      keep2 <- if (twomode) p[lab2] else p[lab1]
      if (nrow(comp) == 0 || isTRUE(cs$nocomp)) keep2 <- rep(TRUE, length(lab2))
      if (!is.null(oppmat)) keep2 <- keep2 & oppmat[k, ]
      cand <- if (twomode) keep2 else keep2 & (lab2 != lab1[s])
      if (length(v) == length(lab2)) L <- lab2 else L <- lab2[cand]
      ok <- if (length(v) == length(lab2)) cand else rep(TRUE, length(L))
      rows[[length(rows) + 1]] <- data.frame(event = k, actor = L[ok], value = v[ok], stringsAsFactors = FALSE)
    } else {
      keep1 <- if (nrow(comp) == 0 || isTRUE(cs$nocomp)) rep(TRUE, length(lab1)) else p[lab1]
      L <- if (length(v) == length(lab1)) lab1 else lab1[keep1]
      ok <- if (length(v) == length(lab1)) keep1 else rep(TRUE, length(L))
      rows[[length(rows) + 1]] <- data.frame(event = k, actor = L[ok], value = v[ok], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# outliers (goldfish's examineOutliers; the settings), change points (the package changepoint as examineChangepoints calls it)
refseries <- function(f, id, cdf_sub) {
  m <- f$m; prep <- f$prep
  x <- m$intervalLogL[prep$orderEvents == 1]
  N <- length(x)
  tm <- prep$eventTime[prep$orderEvents == 1]
  ol <- list(); cl <- list()
  sets <- list(list("Top", 3, NA), list("Top", 1, NA), list("IQR", 3, NA), list("IQR", 1, NA), list("Hampel", 3, NA), list("Hampel", 3, max(1, floor(N / 5))),
               list("Hampel", 2, 3), list("Hampel", 4, max(2, floor(N / 3))))
  for (s in sets) {
    r <- tryCatch({
      p <- if (is.na(s[[3]])) examineOutliers(m, method = s[[1]], parameter = s[[2]]) else examineOutliers(m, method = s[[1]], parameter = s[[2]], window = s[[3]])
      if (inherits(p, "ggplot")) paste(which(p$data$outlier == "YES"), collapse = " ") else ""
    }, error = function(e) "ERR")
    ol[[length(ol) + 1]] <- data.frame(case = id, method = s[[1]], parameter = s[[2]], window = ifelse(is.na(s[[3]]), -1, s[[3]]), idx = r, stringsAsFactors = FALSE)
  }
  win0 <- max(table(tm))
  for (mom in c("mean", "variance")) for (meth in c("PELT", "AMOC", "BinSeg")) for (w in c(NA, 3, 5)) {
    wd <- if (is.na(w)) win0 else w
    r <- tryCatch({
      cpt <- if (mom == "mean") changepoint::cpt.mean(x, method = meth, minseglen = wd) else changepoint::cpt.var(x, method = meth, minseglen = wd)
      pts <- attributes(cpt)$cpts
      if (anyDuplicated(tm[pts])) pts <- pts[!duplicated(tm[pts], fromLast = TRUE)]
      paste(pts, collapse = " ")
    }, error = function(e) "ERR", warning = function(w) { cpt <- NULL; "WARN" })
    if (identical(r, "WARN")) {
      r <- tryCatch(suppressWarnings({
        cpt <- if (mom == "mean") changepoint::cpt.mean(x, method = meth, minseglen = wd) else changepoint::cpt.var(x, method = meth, minseglen = wd)
        pts <- attributes(cpt)$cpts
        if (anyDuplicated(tm[pts])) pts <- pts[!duplicated(tm[pts], fromLast = TRUE)]
        paste(pts, collapse = " ")
      }), error = function(e) "ERR")
    }
    cl[[length(cl) + 1]] <- data.frame(case = id, moment = mom, method = meth, window = ifelse(is.na(w), -1, w), idx = r, stringsAsFactors = FALSE)
  }
  # examineChangepoints itself must give the same change points as the direct call (checked for the default window)
  chk <- tryCatch({
    p <- examineChangepoints(m, moment = "mean", method = "PELT")
    if (inherits(p, "ggplot")) {
      v <- p$layers[[3]]$data$xintercept; if (is.null(v)) v <- p$layers[[3]]$geom_params$xintercept
      v
    } else NULL
  }, error = function(e) NA)
  list(ol = do.call(rbind, ol), cl = do.call(rbind, cl), chk = chk, tm = tm)
}

res <- NULL; cdf <- NULL; icdf <- NULL; oldf <- NULL; cpdf <- NULL
record <- function(id, f, cs, sub, opts, directed, twomode = FALSE, oppmat = NULL, tolb = 1e-6, tolse = 1e-6, exo = "") {
  m <- f$m; prep <- f$prep
  ll <- m$intervalLogL
  write.csv(data.frame(type = prep$orderEvents, ll = ll, time = c(prep$eventTime)), file.path(out, paste0("ref_ll_", id, ".csv")), row.names = FALSE)
  pr <- refprobs(f, cs, sub, twomode, oppmat)
  write.csv(pr, file.path(out, paste0("ref_pr_", id, ".csv")), row.names = FALSE)
  icdf <<- rbind(icdf, data.frame(case = id, ll = m$logLikelihood, nEvents = m$nEvents, nParams = m$nParams, aic = AIC(m), bic = BIC(m),
                                  avg = as.numeric(logLik(m, avgPerEvent = TRUE)), nDep = sum(prep$orderEvents == 1)))
  s <- refseries(f, id, sub)
  oldf <<- rbind(oldf, s$ol); cpdf <<- rbind(cpdf, s$cl)
  if (!is.null(s$chk) && !identical(s$chk, NA)) {
    dflt <- s$cl[s$cl$moment == "mean" & s$cl$method == "PELT" & s$cl$window == -1, "idx"]
    pts <- as.numeric(strsplit(dflt, " ")[[1]])
    if (length(pts) > 1 || (length(pts) == 1 && s$tm[pts] != max(s$tm))) {
      if (!isTRUE(all.equal(sort(unique(s$tm[pts])), sort(unique(as.numeric(s$chk)))))) cat(id, "WARNING: examineChangepoints differs from the direct call\n")
    }
  }
  res <<- rbind(res, data.frame(case = id, term = c(names(f$b), "ll"), estimate = c(f$b, f$ll), se = c(f$se, NA)))
  cdf <<- rbind(cdf, data.frame(id = id, sub = sub, opts = opts, directed = directed, tolb = tolb, tolse = tolse, exo = exo))
  cat(id, sub, ":", paste(sprintf("%.6f", c(f$b, f$ll)), collapse = " "), " nEvents", m$nEvents, "dep", sum(prep$orderEvents == 1), "\n")
}
stable <- function(f, sub) {
  !(is.null(f) || any(!is.finite(f$se)) || max(abs(f$b[names(f$b) != "Intercept"])) > 8 || max(f$se) > (if (sub == "choice_coordination") 0.6 else 3))
}
writecase <- function(cs, id, twomode = FALSE) {
  write.csv(cs$ev[, c("time", "sender", "receiver", "increment")], file.path(out, paste0("events_", id, ".csv")), row.names = FALSE)
  if (twomode) nd <- data.frame(label = cs$lab, x = 0, y = 0, g = 1, present0 = as.integer(cs$init))
  else { nd <- cs$nodes[, c("label", "x", "y", "g")]; nd$present0 <- as.integer(cs$nodes$present) }
  write.csv(nd, file.path(out, paste0("nodes_", id, ".csv")), row.names = FALSE)
  if (!isTRUE(cs$nocomp)) write.csv(data.frame(time = cs$comp$time, node = cs$comp$node, value = as.integer(cs$comp$value)), file.path(out, paste0("comp_", id, ".csv")), row.names = FALSE)
  if (!is.null(cs$chx)) write.csv(cs$chx, file.path(out, paste0("changes_x_", id, ".csv")), row.names = FALSE)
}

# ---------------------------------------------------------------- family r: random combinations
P <- function(f, o, kind = "", v = "") list(f = f, o = o, kind = kind, v = v)
pieces <- list(
  choice = list(P("inertia(net%W)", "inertia", "iw"), P("recip(net%W)", "recip", "rw"), P("indeg(net)", "indeg"),
    P("same(nodes$%s)", "same(%s)", "v", "g"), P("diff(nodes$%s)", "diff(%s)", "v", "xy"), P("sim(nodes$%s)", "sim(%s)", "v", "xy"),
    P("alter(nodes$%s)", "alter(%s)", "v", "xy"), P("outdeg(net)", "outdeg"),
    P("trans(net)", "trans"), P("nodeTrans(net)", "nodetrans"), P("egoAlterInt(nodes$x, nodes$y)", "egoalterint(x y)", "xy"), P("tertius(net, nodes$%s)", "tertius(%s)", "v", "xy")),
  rate = list(P("indeg(net%W)", "indeg", "xw"), P("outdeg(net%W)", "outdeg", "ow"), P("ego(nodes$%s)", "ego(%s)", "v", "xy"),
    P("tertius(net, nodes$%s)", "tertius(%s)", "v", "xy")),
  choice_coordination = list(P("inertia(net%W)", "inertia", "iw"), P("indeg(net)", "indeg"), P("same(nodes$%s)", "same(%s)", "v", "g"),
    P("diff(nodes$%s)", "diff(%s)", "v", "xy"), P("sim(nodes$%s)", "sim(%s)", "v", "xy"), P("alter(nodes$%s)", "alter(%s)", "v", "xy"),
    P("tertius(net, nodes$%s)", "tertius(%s)", "v", "xy"), P("egoAlterInt(nodes$x, nodes$y)", "egoalterint(x y)", "xy")))
mkspec <- function(sub, icpt, p_win, p_tvc) {
  pc <- pieces[[sub]]
  take <- which(runif(length(pc)) < 0.4)
  if (length(take) == 0) take <- sample(length(pc), 1)
  take <- take[seq_len(min(length(take), 4))]
  forms <- character(); opts <- character(); used <- character()
  for (t in take) {
    pp <- pc[[t]]
    if (pp$kind %in% c("iw", "rw", "xw", "ow")) {
      w <- if (runif(1) < p_win) sample(c(4, 6, 9, 14), 1) else NA
      forms <- c(forms, sub("%W", if (is.na(w)) "" else paste0(", window = ", w), pp$f, fixed = TRUE))
      o <- pp$o
      if (!is.na(w)) o <- paste0(o, " ", switch(pp$kind, iw = "inertiawindow", rw = "recipwindow", xw = "indegwindow", ow = "outdegwindow"), "(", w, ")")
      opts <- c(opts, o)
    } else if (pp$kind == "v") {
      v <- sample(strsplit(pp$v, "")[[1]], 1)
      forms <- c(forms, sprintf(sub("%s", "%s", pp$f, fixed = TRUE), v)); opts <- c(opts, sprintf(pp$o, v)); used <- c(used, v)
    } else if (pp$kind == "xy") { forms <- c(forms, pp$f); opts <- c(opts, pp$o); used <- c(used, "x", "y") }
    else { forms <- c(forms, pp$f); opts <- c(opts, pp$o) }
  }
  tvc <- ("x" %in% used) && runif(1) < p_tvc
  list(form = paste(c(if (icpt) "1", forms), collapse = " + "), opts = paste(c(opts, if (icpt) "intercept", if (tvc) "change(x chg_x)"), collapse = " "), tvc = tvc)
}
set.seed(2027)
nrand <- 72
for (k in 1:nrand) {
  sub <- c("choice", "rate", "choice", "rate", "choice_coordination", "choice")[(k - 1) %% 6 + 1]
  icpt <- (sub == "rate" && k %% 4 == 2)
  sp <- mkspec(sub, icpt, 0.4, 0.5)
  nocomp <- (k %% 3 == 0)
  id <- sprintf("r%02d", k)
  cs <- mkcase(900 + k, sample(6:12, 1), sample(40:110, 1), pabs = sample(c(0.1, 0.25, 0.45), 1), ncomp = sample(6:30, 1), tvc = sp$tvc,
               directed = (sub != "choice_coordination"), nocomp = nocomp)
  f <- tryCatch(fitgf(cs, sp$form, sub), error = function(e) { cat(id, "SKIPPED", conditionMessage(e), "\n"); NULL })
  if (!stable(f, sub)) { if (!is.null(f)) cat(id, "SKIPPED (unstable)\n"); next }
  writecase(cs, id)
  record(id, f, cs, sub, paste(sp$opts, if (!nocomp) "present(present0) compchange(comp)" else ""), as.integer(sub != "choice_coordination"),
         tolb = if (sub == "choice_coordination") 1e-5 else 1e-6, tolse = if (sub == "choice_coordination") 0.1 else 1e-6)
}

# ---------------------------------------------------------------- family m: two-mode
two <- list(
  list(id = "m01", seed = 61, n1 = 6, n2 = 5, N = 60, form = "inertia(net) + indeg(net)", opts = "inertia indeg", sub = "choice"),
  list(id = "m02", seed = 62, n1 = 7, n2 = 5, N = 60, form = "outdeg(net)", opts = "outdeg", sub = "rate"),
  list(id = "m03", seed = 63, n1 = 7, n2 = 5, N = 60, form = "1 + outdeg(net)", opts = "outdeg intercept", sub = "rate"),
  list(id = "m04", seed = 64, n1 = 8, n2 = 6, N = 80, pabs = 0.4, ncomp = 20, form = "inertia(net)", opts = "inertia", sub = "choice"),
  list(id = "m05", seed = 65, n1 = 8, n2 = 6, N = 80, pabs = 0.4, ncomp = 20, form = "1 + outdeg(net)", opts = "outdeg intercept", sub = "rate"),
  list(id = "m06", seed = 66, n1 = 8, n2 = 6, N = 80, pabs = 0.1, ncomp = 20, form = "indeg(net)", opts = "indeg", sub = "choice"))
for (cc in two) {
  cs <- mk2(cc$seed, cc$n1, cc$n2, cc$N, pabs = if (is.null(cc$pabs)) 0.25 else cc$pabs, ncomp = if (is.null(cc$ncomp)) 10 else cc$ncomp)
  f <- tryCatch(fit2(cs, cc$form, cc$sub), error = function(e) { cat(cc$id, "SKIPPED", conditionMessage(e), "\n"); NULL })
  if (!stable(f, cc$sub)) { if (!is.null(f)) cat(cc$id, "SKIPPED (unstable)\n"); next }
  cs$nocomp <- FALSE
  writecase(cs, cc$id, twomode = TRUE)
  record(cc$id, f, cs, cc$sub, paste(cc$opts, "present(present0) compchange(comp)"), 2L, twomode = TRUE)
}

# ---------------------------------------------------------------- family o: choice with an opportunity list per event
for (k in 1:6) {
  id <- sprintf("o%02d", k)
  set.seed(5000 + k)
  n <- sample(6:9, 1)
  cs <- mkcase(5000 + k, n, sample(35:70, 1), tvc = (k %% 2 == 0), nocomp = TRUE, ties = (k %% 3 != 0), cycle = TRUE)
  # event times are the sequence numbers (opportunities refer to the event number; nwdynam's event order is the chronological order)
  forms <- c("inertia(net) + indeg(net)", "inertia(net) + alter(nodes$x)", "recip(net) + indeg(net) + sim(nodes$x)", "inertia(net) + same(nodes$g)", "indeg(net) + outdeg(net)", "inertia(net) + recip(net)")[k]
  opts <- c("inertia indeg", "inertia alter(x)", "recip indeg sim(x)", "inertia same(g)", "indeg outdeg", "inertia recip")[k]
  if (k %% 2 == 0 && !grepl("x", forms)) { cs$chx <- NULL }
  if (k %% 2 == 0 && grepl("x", forms)) { } else cs$chx <- NULL
  oppl <- lapply(1:cs$N, function(i) {
    s <- match(cs$ev$sender[i], cs$lab); r <- match(cs$ev$receiver[i], cs$lab)
    pool <- setdiff(1:n, c(s, r))
    nex <- sample(0:min(2, length(pool) - 2), 1)
    ex <- if (nex > 0) sample(pool, nex) else integer(0)
    setdiff(1:n, ex)
  })
  f <- tryCatch(fitgf(cs, forms, "choice", opp = oppl), error = function(e) { cat(id, "SKIPPED", conditionMessage(e), "\n"); NULL })
  if (!stable(f, "choice")) { if (!is.null(f)) cat(id, "SKIPPED (unstable)\n"); next }
  writecase(cs, id)
  av <- do.call(rbind, lapply(seq_along(oppl), function(i) data.frame(event = i, actor = oppl[[i]])))
  write.csv(av, file.path(out, paste0("avail_", id, ".csv")), row.names = FALSE)
  oppmat <- t(vapply(oppl, function(o) (1:n) %in% o, logical(n)))
  record(id, f, cs, "choice", paste(opts, if (!is.null(cs$chx)) "change(x chg_x)" else "", "oppframe(oppdata) opportunities(event actor)"), 1L, oppmat = oppmat)
}

# ---------------------------------------------------------------- family x: netterms (effects on other networks)
mkexo <- function(seed, n = 8, N = 60, exo = list(), kind = "one") {
  set.seed(seed)
  lab <- paste0("A", 1:n)
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) { s <- sample(lab, 1); snd[k] <- s; rcv[k] <- sample(setdiff(lab, s), 1) }
  snd[1:n] <- sample(lab); for (k in 1:n) rcv[k] <- sample(setdiff(lab, snd[k]), 1)
  for (k in seq(5, N, by = 4)) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }
  tm <- round(cumsum(runif(N, 0.3, 3)) * 4) / 4
  for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  tm <- sort(tm)
  ev <- data.frame(time = tm, sender = snd, receiver = rcv, increment = 1, stringsAsFactors = FALSE)
  nodes <- data.frame(label = lab, present = TRUE, x = round(rnorm(n), 2), y = round(runif(n, 0, 2), 2), g = sample(1:3, n, replace = TRUE), stringsAsFactors = FALSE)
  ex <- list()
  for (nm in names(exo)) {
    sp <- exo[[nm]]
    if (sp$type == "stat") {
      M <- matrix(0, n, n)
      for (i in 1:n) for (j in 1:n) if (i != j && runif(1) < 0.3) M[i, j] <- 1
      ex[[nm]] <- list(type = "stat", directed = TRUE, M = M)
    } else {
      ne <- 45
      tt <- sort(round(runif(ne, min(tm) - 2, max(tm) + 1) * 4) / 4)
      tt[sample(ne, floor(ne / 4))] <- sample(tm, floor(ne / 4), TRUE)
      tt <- sort(tt)
      cur <- matrix(0, n, n)
      s2 <- character(ne); r2 <- character(ne); v2 <- numeric(ne)
      for (k in 1:ne) {
        i <- sample(n, 1); j <- sample(setdiff(1:n, i), 1)
        if (cur[i, j] > 0 && runif(1) < 0.3) v <- -cur[i, j] else v <- 1
        cur[i, j] <- cur[i, j] + v
        s2[k] <- lab[i]; r2[k] <- lab[j]; v2[k] <- v
      }
      ex[[nm]] <- list(type = "dyn", directed = TRUE, ev = data.frame(time = tt, sender = s2, receiver = r2, increment = v2, stringsAsFactors = FALSE))
    }
  }
  list(lab = lab, ev = ev, nodes = nodes, exo = ex, directed = TRUE, comp = data.frame(time = numeric(0), node = character(0), value = logical(0)), chx = NULL, nocomp = TRUE, n = n, N = N)
}
EXO <- list(x1 = list(type = "dyn"), x2 = list(type = "stat"))
xspec <- list(
  list("choice", "tie(x1)", "netterms(tie(x1))", c("x1")), list("choice", "tie(x2)", "netterms(tie(x2))", "x2"),
  list("choice", "inertia(net) + trans(x1)", "inertia netterms(trans(x1))", "x1"), list("choice", "indeg(net) + cycle(x2)", "indeg netterms(cycle(x2))", "x2"),
  list("choice", "inertia(net) + mixedTrans(list(net, x1), netUpdate = 1)", "inertia netterms(mixedtrans(net x1))", "x1"),
  list("choice", "tie(x1, window = 6) + four(x2)", "netterms(tie(x1, window(6)) four(x2))", c("x1", "x2")),
  list("choice", "inertia(net) + commonSender(x1)", "inertia netterms(commonsender(x1))", "x1"),
  list("choice_coordination", "inertia(net) + tie(x2)", "inertia netterms(tie(x2))", "x2"),
  list("choice_coordination", "indeg(net) + tie(x1)", "indeg netterms(tie(x1))", "x1"),
  list("choice", "inertia(net) + tie(x1) + trans(x2)", "inertia netterms(tie(x1) trans(x2))", c("x1", "x2")))
for (k in seq_along(xspec)) {
  sp <- xspec[[k]]; id <- sprintf("x%02d", k)
  cs <- mkexo(7000 + k, n = sample(6:8, 1), N = sample(45:70, 1), exo = EXO[sp[[4]]])
  coord <- (sp[[1]] == "choice_coordination")
  cs$directed <- !coord
  if (coord) for (nm in names(cs$exo)) if (cs$exo[[nm]]$type == "stat") { M <- cs$exo[[nm]]$M; M[lower.tri(M)] <- t(M)[lower.tri(M)]; cs$exo[[nm]]$M <- M; cs$exo[[nm]]$directed <- FALSE } else cs$exo[[nm]]$directed <- FALSE
  f <- tryCatch(fitgf(cs, sp[[2]], sp[[1]]), error = function(e) { cat(id, "SKIPPED", conditionMessage(e), "\n"); NULL })
  if (!stable(f, sp[[1]])) { if (!is.null(f)) cat(id, "SKIPPED (unstable)\n"); next }
  writecase(cs, id)
  exospec <- character()
  for (nm in names(cs$exo)) {
    e <- cs$exo[[nm]]
    if (e$type == "stat") { ij <- which(e$M != 0, arr.ind = TRUE); write.csv(data.frame(sender = ij[, 1], receiver = ij[, 2], value = e$M[ij]), file.path(out, paste0("exo_", id, "_", nm, ".csv")), row.names = FALSE) }
    else write.csv(e$ev, file.path(out, paste0("exo_", id, "_", nm, ".csv")), row.names = FALSE)
    exospec <- c(exospec, paste0(nm, "=", e$type, if (e$directed) "dir" else "und"))
  }
  record(id, f, cs, sp[[1]], sp[[3]], as.integer(!coord), tolb = if (coord) 5e-4 else 1e-6, tolse = if (coord) 0.1 else 1e-6, exo = paste(exospec, collapse = " "))
}

if (!"exo" %in% names(cdf)) cdf$exo <- ""
cdf$exo[is.na(cdf$exo)] <- ""
write.csv(res, file.path(out, "results.csv"), row.names = FALSE)
write.csv(cdf, file.path(out, "cases.csv"), row.names = FALSE)
write.csv(icdf, file.path(out, "post_ic.csv"), row.names = FALSE)
write.csv(oldf, file.path(out, "post_outliers.csv"), row.names = FALSE)
write.csv(cpdf, file.path(out, "post_cpts.csv"), row.names = FALSE)
cat("cases:", nrow(cdf), "\n")
