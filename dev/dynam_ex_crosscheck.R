# nwdynam: effects on other networks, mixed closures, tertiusDiff, tie() of a changing network, ignoreRep against goldfish (R 4.x, goldfish 1.6.12).
# Usage: Rscript dev/dynam_ex_crosscheck.R <output directory> [sections, e.g. a b c d e f g]
# Sections: a = tie() of exogenous networks (static, changing, weighted, windowed, ignoreRep), b = trans/cycle/commonSender/commonReceiver
# through exogenous networks, c = mixed closures, d = tertiusDiff, e = ignoreRep on the effects of the dependent network,
# f = choice_coordination (tie, trans, mixedTrans, tertiusDiff, ignoreRep), g = hand-made edge cases.
# For every case it writes events_<case>.csv (time, sender, receiver, increment), nodes_<case>.csv (label, x, y, g), exo_<case>_<name>.csv (the
# exogenous networks: changing = time, sender, receiver, increment; static = sender, receiver, value), results_<sec>.csv (case, term, estimate, se;
# the last row of a case is "ll") and cases_<sec>.csv (case, sub, options of nwdynam, kind, exogenous networks, tolerances).
# Event times are multiples of 0.25 (exact in the float times of the core); actors are A1..A9 (text order = numeric order).
suppressMessages(library(goldfish))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "."
secs <- if (length(args) >= 2) args[-1] else c("a", "b", "c", "d", "e", "f", "g")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- one random case. exo: named list; each element list(type = "dyn" | "stat", weighted, directed, nev, density)
mkcase <- function(seed, n = 8, N = 60, kind = "one", weights = FALSE, inttime = FALSE, ties = TRUE, exo = list(), canon = FALSE) {
  set.seed(seed)
  lab <- paste0("A", 1:n)
  snd <- character(N); rcv <- character(N)
  for (k in 1:N) { s <- sample(lab, 1); snd[k] <- s; rcv[k] <- sample(setdiff(lab, s), 1) }
  snd[1:n] <- sample(lab); for (k in 1:n) rcv[k] <- sample(setdiff(lab, snd[k]), 1)   # every actor has an event (the network has n actors)
  for (k in seq(5, N, by = 4)) { snd[k] <- snd[k - 3]; rcv[k] <- rcv[k - 3] }
  if (inttime) { tm <- cumsum(sample(0:3, N, TRUE, prob = c(.15, .35, .3, .2))); tm <- tm + 1 }
  else {
    tm <- round(cumsum(runif(N, 0.3, 3)) * 4) / 4
    if (ties) for (k in seq(7, N, by = 9)) tm[k] <- tm[k - 1]
  }
  tm <- sort(tm)
  val <- if (weights) sample(1:3, N, TRUE, prob = c(.6, .25, .15)) else rep(1, N)
  if (kind == "und" && canon) { sw <- match(snd, lab) > match(rcv, lab); tmp <- snd[sw]; snd[sw] <- rcv[sw]; rcv[sw] <- tmp }
  ev <- data.frame(time = tm, sender = snd, receiver = rcv, increment = val, stringsAsFactors = FALSE)
  nodes <- data.frame(label = lab, present = TRUE, x = round(rnorm(n), 2), y = round(runif(n, 0, 2), 2), g = sample(1:3, n, replace = TRUE), stringsAsFactors = FALSE)
  ex <- list()
  for (nm in names(exo)) {
    sp <- exo[[nm]]
    dir <- if (is.null(sp$directed)) TRUE else sp$directed
    wt <- isTRUE(sp$weighted)
    if (sp$type == "stat") {
      dens <- if (is.null(sp$density)) 0.3 else sp$density
      M <- matrix(0, n, n)
      for (i in 1:n) for (j in 1:n) if (i != j && runif(1) < dens) M[i, j] <- if (wt) sample(1:3, 1) else 1
      if (!dir) { M[lower.tri(M)] <- t(M)[lower.tri(M)] }
      ex[[nm]] <- list(type = "stat", directed = dir, M = M)
    } else {
      ne <- if (is.null(sp$nev)) 40 else sp$nev
      tt <- sort(round(runif(ne, min(tm) - 2, max(tm) + 1) * 4) / 4)
      if (ne >= 6) tt[sample(ne, floor(ne / 4))] <- sample(tm, floor(ne / 4), TRUE)   # some at the times of dependent events
      tt <- sort(tt)
      cur <- matrix(0, n, n)
      s2 <- character(ne); r2 <- character(ne); v2 <- numeric(ne)
      for (k in 1:ne) {
        i <- sample(n, 1); j <- sample(setdiff(1:n, i), 1)
        if (cur[i, j] > 0 && runif(1) < 0.3) v <- -cur[i, j] else v <- if (wt) sample(1:3, 1) else (if (cur[i, j] > 0) 1 else 1)
        cur[i, j] <- cur[i, j] + v; if (!dir) cur[j, i] <- cur[j, i] + v
        s2[k] <- lab[i]; r2[k] <- lab[j]; v2[k] <- v
      }
      ex[[nm]] <- list(type = "dyn", directed = dir, ev = data.frame(time = tt, sender = s2, receiver = r2, increment = v2, stringsAsFactors = FALSE))
    }
  }
  list(lab = lab, kind = kind, ev = ev, nodes = nodes, exo = ex, n = n)
}

# ---- fit in goldfish. form: formula with objects net, nd and the exogenous networks by name
fitgf <- function(cs, form, sub) {
  ev <- cs$ev
  assign("ev", ev, globalenv())
  nd <- defineNodes(cs$nodes)
  assign("nd", nd, globalenv())
  net <- defineNetwork(nodes = nd, directed = (cs$kind == "one"))
  net <- linkEvents(net, changeEvent = ev, nodes = nd)
  dep <- defineDependentEvents(events = ev, nodes = nd, defaultNetwork = net)
  assign("net", net, globalenv()); assign("dep", dep, globalenv())
  for (nm in names(cs$exo)) {
    e <- cs$exo[[nm]]
    if (e$type == "stat") {
      x <- defineNetwork(matrix = e$M, nodes = nd, directed = e$directed)
    } else {
      evn <- paste0("exoev_", nm); assign(evn, e$ev, globalenv())
      x <- defineNetwork(nodes = nd, directed = e$directed)
      x <- eval(substitute(linkEvents(x, changeEvent = EV, nodes = nd), list(EV = as.name(evn))))
    }
    assign(nm, x, globalenv())
  }
  ei <- list(maxIterations = 200, initialDamping = 10, maxScoreStopCriterion = 1e-9)
  m <- suppressWarnings(estimate(as.formula(paste("dep ~", form), env = globalenv()), model = "DyNAM", subModel = sub, estimationInit = ei))
  list(b = coef(m), se = m$standardErrors, ll = as.numeric(m$logLikelihood))
}

res <- list(); cdf <- list()
runspec <- function(sp) {
  cs <- do.call(mkcase, c(list(seed = sp$seed), sp$mk))
  cs$exo <- cs$exo[vapply(names(cs$exo), function(nm) grepl(paste0("\\b", nm, "\\b"), sp$form), logical(1))]
  f <- tryCatch(fitgf(cs, sp$form, sp$sub), error = function(e) { cat(sp$id, "SKIPPED:", sp$form, "|", conditionMessage(e), "\n"); NULL })
  if (is.null(f)) return(invisible())
  if (any(!is.finite(f$se)) || max(abs(f$b)) > 8 || max(f$se) > (if (sp$sub == "choice_coordination") 0.6 else 3)) { cat(sp$id, "SKIPPED (unstable)\n"); return(invisible()) }
  id <- sp$id
  write.csv(cs$ev[, c("time", "sender", "receiver", "increment")], file.path(out, paste0("events_", id, ".csv")), row.names = FALSE)
  write.csv(cs$nodes[, c("label", "x", "y", "g")], file.path(out, paste0("nodes_", id, ".csv")), row.names = FALSE)
  exospec <- character()
  for (nm in names(cs$exo)) {
    e <- cs$exo[[nm]]
    if (e$type == "stat") {
      ij <- which(e$M != 0, arr.ind = TRUE)
      write.csv(data.frame(sender = ij[, 1], receiver = ij[, 2], value = e$M[ij]), file.path(out, paste0("exo_", id, "_", nm, ".csv")), row.names = FALSE)
    } else write.csv(e$ev, file.path(out, paste0("exo_", id, "_", nm, ".csv")), row.names = FALSE)
    exospec <- c(exospec, paste0(nm, "=", e$type, if (e$directed) "dir" else "und"))
  }
  cat(id, ":", paste(sprintf("%.8f", c(f$b, f$ll)), collapse = " "), "\n")
  coord <- (sp$sub == "choice_coordination")
  res[[length(res) + 1]] <<- data.frame(sec = sp$sec, case = id, term = c(names(f$b), "ll"), estimate = c(f$b, f$ll), se = c(f$se, NA))
  cdf[[length(cdf) + 1]] <<- data.frame(sec = sp$sec, id = id, sub = sp$sub, opts = sp$opts, kind = cs$kind, exo = paste(exospec, collapse = " "),
    tolb = if (coord) 5e-4 else 1e-6, tolse = if (coord) 0.1 else 1e-6)
}

# pieces: list(goldfish term, nwdynam text, kind) with kind "b" = option of nwdynam, "x" = term of netterms(); placeholders %W (window), %s (x or y)
randspecs <- function(sec, pieces, nrand, seed0, sub, mkargs, maxeff = 3, idoff = 0, must = NULL, groups = NULL) {
  set.seed(seed0)
  sp <- list()
  for (k in 1:nrand) {
    take <- which(runif(length(pieces)) < 0.5)
    if (length(take) == 0) take <- sample(length(pieces), 1)
    if (!is.null(groups)) take <- take[!duplicated(groups[take])]
    take <- take[seq_len(min(length(take), maxeff))]
    forms <- character(); bopts <- character(); xopts <- character()
    for (t in take) {
      pp <- pieces[[t]]
      v <- if (length(pp) >= 4) sample(strsplit(pp[[4]], "")[[1]], 1) else ""
      f <- gsub("%s", v, pp[[1]], fixed = TRUE); o <- gsub("%s", v, pp[[2]], fixed = TRUE)
      forms <- c(forms, f)
      if (pp[[3]] == "b") bopts <- c(bopts, o) else xopts <- c(xopts, o)
    }
    if (!is.null(must)) { forms <- c(forms, must[[1]]); if (must[[3]] == "b") bopts <- c(bopts, must[[2]]) else xopts <- c(xopts, must[[2]]) }
    # nwdynam's order: the options first, then the netterms(); the formula must follow it
    isb <- vapply(c(take, if (!is.null(must)) 0L else integer()), function(t) { if (t == 0) must[[3]] == "b" else pieces[[t]][[3]] == "b" }, logical(1))
    forms <- c(forms[isb], forms[!isb])
    fm <- paste(forms, collapse = " + ")
    op <- paste(c(bopts, if (length(xopts)) paste0("netterms(", paste(xopts, collapse = " "), ")")), collapse = " ")
    while (grepl("%W", fm, fixed = TRUE)) { w <- sample(c(2, 3, 4.5, 6, 9), 1); fm <- sub("%W", w, fm, fixed = TRUE); op <- sub("%W", w, op, fixed = TRUE) }
    sp[[k]] <- list(sec = sec, id = sprintf("%s%02d", sec, k + idoff), seed = seed0 * 100 + k, sub = sub, mk = mkargs(k), form = fm, opts = op)
  }
  sp
}
one <- function(sec, id, seed, sub, mk, form, opts) list(sec = sec, id = id, seed = seed, sub = sub, mk = mk, form = form, opts = opts)
EXO <- list(x1 = list(type = "dyn", directed = TRUE, nev = 45), x2 = list(type = "stat", directed = TRUE), x3 = list(type = "dyn", directed = TRUE, weighted = TRUE, nev = 45),
            x4 = list(type = "stat", directed = TRUE, weighted = TRUE))
exo_for <- function(k, names) EXO[names]

specs <- list()
# ================= a: tie() of exogenous networks =================
if ("a" %in% secs) {
  pa <- list(list("tie(x1)", "tie(x1)", "x"), list("tie(x2)", "tie(x2)", "x"), list("tie(x3, weighted = TRUE)", "tie(x3, weighted)", "x"),
             list("tie(x4, weighted = TRUE)", "tie(x4, weighted)", "x"), list("tie(x1, window = %W)", "tie(x1, window(%W))", "x"),
             list("tie(x3, weighted = TRUE, window = %W)", "tie(x3, weighted window(%W))", "x"), list("tie(x1, ignoreRep = TRUE)", "tie(x1, ignorerep)", "x"))
  ba <- list(list("inertia(net)", "inertia", "b"), list("indeg(net)", "indeg", "b"))
  sa <- list(one("a", "a01", 1, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "tie(x1)", "netterms(tie(x1))"),
             one("a", "a02", 2, "choice", list(n = 7, N = 50, exo = EXO["x3"]), "tie(x3, weighted = TRUE)", "netterms(tie(x3, weighted))"),
             one("a", "a03", 3, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "tie(x1, window = 6)", "netterms(tie(x1, window(6)))"),
             one("a", "a04", 4, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "tie(x2)", "netterms(tie(x2))"),
             one("a", "a05", 5, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "inertia(net) + tie(x1, ignoreRep = TRUE)", "inertia netterms(tie(x1, ignorerep))"))
  sa2 <- randspecs("a", c(ba, pa), 40, 11, "choice", function(k) list(n = sample(6:9, 1), N = sample(40:80, 1), exo = EXO, weights = (k %% 4 == 0),
         inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 10)
  specs <- c(specs, sa, sa2)
}
# ================= b: closure effects through exogenous networks =================
if ("b" %in% secs) {
  pb <- list(list("trans(x1)", "trans(x1)", "x"), list("trans(x2)", "trans(x2)", "x"), list("cycle(x1)", "cycle(x1)", "x"), list("cycle(x2)", "cycle(x2)", "x"),
             list("commonSender(x1)", "commonsender(x1)", "x"), list("commonSender(x2)", "commonsender(x2)", "x"),
             list("commonReceiver(x1)", "commonreceiver(x1)", "x"), list("commonReceiver(x2)", "commonreceiver(x2)", "x"),
             list("trans(x1, window = %W)", "trans(x1, window(%W))", "x"), list("cycle(net, window = %W)", "cycle(net, window(%W))", "x"),
             list("trans(net, ignoreRep = TRUE)", "trans(net, ignorerep)", "x"), list("commonSender(x1, ignoreRep = TRUE)", "commonsender(x1, ignorerep)", "x"),
             list("four(x1)", "four(x1)", "x"), list("four(x2)", "four(x2)", "x"), list("four(x1, window = %W)", "four(x1, window(%W))", "x"), list("four(net, window = %W)", "four(net, window(%W))", "x"),
             list("four(x1, ignoreRep = TRUE)", "four(x1, ignorerep)", "x"))
  bb <- list(list("inertia(net)", "inertia", "b"), list("indeg(net)", "indeg", "b"))
  sb <- list(one("b", "b01", 1, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "trans(x1)", "netterms(trans(x1))"),
             one("b", "b02", 2, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "cycle(x2)", "netterms(cycle(x2))"),
             one("b", "b03", 3, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "commonSender(x1)", "netterms(commonsender(x1))"),
             one("b", "b04", 4, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "commonReceiver(x2)", "netterms(commonreceiver(x2))"),
             one("b", "b05", 5, "choice", list(n = 7, N = 50), "trans(net, window = 6)", "netterms(trans(net, window(6)))"),
             one("b", "b06", 6, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "four(x1)", "netterms(four(x1))"),
             one("b", "b07", 7, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "four(x2)", "netterms(four(x2))"))
  sb2 <- randspecs("b", c(bb, pb), 45, 21, "choice", function(k) list(n = sample(6:9, 1), N = sample(40:80, 1), exo = EXO, weights = (k %% 4 == 0),
         inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 10)
  specs <- c(specs, sb, sb2)
}
# ================= c: mixed closures =================
if ("c" %in% secs) {
  pc <- list(list("mixedTrans(list(net, x1), netUpdate = 1)", "mixedtrans(net x1)", "x"), list("mixedTrans(list(x1, net), netUpdate = 2)", "mixedtrans(x1 net)", "x"),
             list("mixedTrans(list(net, x2), netUpdate = 1)", "mixedtrans(net x2)", "x"), list("mixedTrans(list(x2, net), netUpdate = 2)", "mixedtrans(x2 net)", "x"),
             list("mixedCycle(list(net, x1), netUpdate = 1)", "mixedcycle(net x1)", "x"), list("mixedCycle(list(x2, net), netUpdate = 2)", "mixedcycle(x2 net)", "x"),
             list("mixedCommonSender(list(net, x1), netUpdate = 1)", "mixedcommonsender(net x1)", "x"), list("mixedCommonSender(list(x1, net), netUpdate = 2)", "mixedcommonsender(x1 net)", "x"),
             list("mixedCommonReceiver(list(net, x2), netUpdate = 1)", "mixedcommonreceiver(net x2)", "x"), list("mixedCommonReceiver(list(x2, net), netUpdate = 2)", "mixedcommonreceiver(x2 net)", "x"),
             list("mixedTrans(list(net, x1), netUpdate = 1, window = %W)", "mixedtrans(net x1, window(%W))", "x"),
             list("mixedTrans(list(net, x2), netUpdate = 1, ignoreRep = TRUE)", "mixedtrans(net x2, ignorerep)", "x"))
  bc <- list(list("inertia(net)", "inertia", "b"), list("indeg(net)", "indeg", "b"))
  sc <- list(one("c", "c01", 1, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "mixedTrans(list(net, x1), netUpdate = 1)", "netterms(mixedtrans(net x1))"),
             one("c", "c02", 2, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "mixedTrans(list(x1, net), netUpdate = 2)", "netterms(mixedtrans(x1 net))"),
             one("c", "c03", 3, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "mixedCommonSender(list(net, x2), netUpdate = 1)", "netterms(mixedcommonsender(net x2))"),
             one("c", "c04", 4, "choice", list(n = 7, N = 50, exo = EXO[c("x1", "x2")]), "mixedTrans(list(x1, x2), netUpdate = 1)", "netterms(mixedtrans(x1 x2))"),
             one("c", "c05", 5, "choice", list(n = 7, N = 50, exo = EXO[c("x1", "x2")]), "mixedTrans(list(x1, x2), netUpdate = 2)", "netterms(mixedtrans(x1 x2))"),
             one("c", "c06", 6, "choice", list(n = 7, N = 50, exo = EXO[c("x2", "x4")]), "mixedTrans(list(x2, x4), netUpdate = 1)", "netterms(mixedtrans(x2 x4))"))
  sc2 <- randspecs("c", c(bc, pc), 45, 31, "choice", function(k) list(n = sample(6:9, 1), N = sample(40:80, 1), exo = EXO, weights = (k %% 4 == 0),
         inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 10)
  specs <- c(specs, sc, sc2)
}
# ================= d: tertiusDiff =================
if ("d" %in% secs) {
  pd <- list(list("tertiusDiff(net, nd$%s)", "tertiusdiff(%s)", "x", "xy"), list("tertiusDiff(x1, nd$%s)", "tertiusdiff(%s x1)", "x", "xy"), list("tertiusDiff(x2, nd$%s)", "tertiusdiff(%s x2)", "x", "xy"),
             list("tertiusDiff(net, nd$x, window = %W)", "tertiusdiff(x, window(%W))", "x"), list("tertiusDiff(x1, nd$y, window = %W)", "tertiusdiff(y x1, window(%W))", "x"),
             list("tertiusDiff(net, nd$y, ignoreRep = TRUE)", "tertiusdiff(y, ignorerep)", "x"))
  bd <- list(list("inertia(net)", "inertia", "b"), list("indeg(net)", "indeg", "b"), list("alter(nd$%s)", "alter(%s)", "b", "xy"))
  sd <- list(one("d", "d01", 1, "choice", list(n = 7, N = 50), "tertiusDiff(net, nd$x)", "netterms(tertiusdiff(x))"),
             one("d", "d02", 2, "choice", list(n = 7, N = 50, exo = EXO["x1"]), "tertiusDiff(x1, nd$x)", "netterms(tertiusdiff(x x1))"),
             one("d", "d03", 3, "choice", list(n = 7, N = 50, exo = EXO["x2"]), "tertiusDiff(x2, nd$y)", "netterms(tertiusdiff(y x2))"),
             one("d", "d04", 4, "choice", list(n = 7, N = 50), "tertiusDiff(net, nd$x, window = 6)", "netterms(tertiusdiff(x, window(6)))"),
             one("d", "d05", 5, "choice", list(n = 8, N = 60), "inertia(net) + tertiusDiff(net, nd$x)", "inertia netterms(tertiusdiff(x))"))
  sd2 <- randspecs("d", c(bd, pd), 50, 41, "choice", function(k) list(n = sample(6:9, 1), N = sample(40:80, 1), exo = EXO, weights = (k %% 4 == 0),
         inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 3, idoff = 10)
  specs <- c(specs, sd, sd2)
}
# ================= e: ignoreRep on the effects of the dependent network =================
if ("e" %in% secs) {
  pe <- list(list("inertia(net)", "inertia", "b"), list("recip(net, ignoreRep = TRUE)", "recip", "b", NULL, "recip"), list("indeg(net, ignoreRep = TRUE)", "indeg", "b", NULL, "indeg"),
             list("same(nd$g)", "same(g)", "b"), list("alter(nd$x)", "alter(x)", "b"),
             list("outdeg(net, ignoreRep = TRUE)", "outdeg", "b", NULL, "outdeg"), list("trans(net, ignoreRep = TRUE)", "trans", "b", NULL, "trans"),
             list("cycle(net, ignoreRep = TRUE)", "cycle", "b", NULL, "cycle"), list("commonSender(net, ignoreRep = TRUE)", "commonsender", "b", NULL, "commonsender"),
             list("commonReceiver(net, ignoreRep = TRUE)", "commonreceiver", "b", NULL, "commonreceiver"), list("four(net, ignoreRep = TRUE)", "four", "b", NULL, "four"),
             list("nodeTrans(net, ignoreRep = TRUE)", "nodetrans", "b", NULL, "nodetrans"), list("tertius(net, nd$x, ignoreRep = TRUE)", "tertius(x)", "b", NULL, "tertius"))
  specs <- c(specs, list(one("e", "e01", 1, "choice", list(n = 7, N = 50), "inertia(net) + indeg(net, ignoreRep = TRUE)", "inertia indeg ignorerep(indeg)"),
                         one("e", "e02", 2, "choice", list(n = 7, N = 50), "inertia(net) + trans(net, ignoreRep = TRUE)", "inertia trans ignorerep(trans)")))
  # (these are built by hand below: the pieces carry their ignorerep word as 5th element)
  set.seed(51)
  for (k in 1:40) {
    take <- which(runif(length(pe)) < 0.4); if (length(take) == 0) take <- sample(length(pe), 1)
    take <- take[seq_len(min(length(take), 3))]
    ign <- Filter(Negate(is.null), lapply(pe[take], function(p) if (length(p) >= 5) p[[5]] else NULL))
    if (length(ign) == 0) next
    ign <- unlist(ign)
    forms <- vapply(pe[take], function(p) p[[1]], ""); opts <- vapply(pe[take], function(p) p[[2]], "")
    # non-ignoreRep pieces use the plain term
    forms <- vapply(pe[take], function(p) if (length(p) >= 5) p[[1]] else p[[1]], "")
    specs[[length(specs) + 1]] <- list(sec = "e", id = sprintf("e%02d", k + 10), seed = 51 * 100 + k, sub = "choice",
      mk = list(n = sample(6:9, 1), N = sample(40:80, 1), weights = (k %% 4 == 0), inttime = (k %% 5 == 0), ties = (k %% 3 != 0)),
      form = paste(forms, collapse = " + "), opts = paste(paste(opts, collapse = " "), " ignorerep(", paste(ign, collapse = " "), ")", sep = ""))
  }
}
# ================= f: choice_coordination =================
if ("f" %in% secs) {
  EXOU <- list(x1 = list(type = "dyn", directed = FALSE, nev = 40), x2 = list(type = "stat", directed = FALSE), x3 = list(type = "dyn", directed = FALSE, weighted = TRUE, nev = 40))
  pf <- list(list("tie(x1)", "tie(x1)", "x"), list("tie(x2)", "tie(x2)", "x"), list("tie(x3, weighted = TRUE)", "tie(x3, weighted)", "x"), list("trans(x1)", "trans(x1)", "x"),
             list("trans(x2)", "trans(x2)", "x"), list("mixedTrans(list(net, x1), netUpdate = 1)", "mixedtrans(net x1)", "x"), list("mixedTrans(list(x2, net), netUpdate = 2)", "mixedtrans(x2 net)", "x"),
             list("tertiusDiff(net, nd$%s)", "tertiusdiff(%s)", "x", "xy"), list("tertiusDiff(x1, nd$x)", "tertiusdiff(x x1)", "x"), list("tie(x1, window = %W)", "tie(x1, window(%W))", "x"),
             list("trans(net, ignoreRep = TRUE)", "trans(net, ignorerep)", "x"),
             list("four(x1)", "four(x1)", "x"), list("four(x2)", "four(x2)", "x"), list("four(net, window = %W)", "four(net, window(%W))", "x"))
  bf <- list(list("inertia(net)", "inertia", "b"), list("indeg(net)", "indeg", "b"))
  sf <- list(one("f", "f01", 1, "choice_coordination", list(kind = "und", n = 7, N = 50, exo = EXOU["x1"]), "tie(x1)", "netterms(tie(x1))"),
             one("f", "f02", 2, "choice_coordination", list(kind = "und", n = 7, N = 50, exo = EXOU["x2"]), "tie(x2)", "netterms(tie(x2))"),
             one("f", "f03", 3, "choice_coordination", list(kind = "und", n = 7, N = 50, exo = EXOU["x1"]), "trans(x1)", "netterms(trans(x1))"),
             one("f", "f04", 4, "choice_coordination", list(kind = "und", n = 7, N = 50, exo = EXOU["x1"]), "mixedTrans(list(net, x1), netUpdate = 1)", "netterms(mixedtrans(net x1))"),
             one("f", "f05", 5, "choice_coordination", list(kind = "und", n = 7, N = 50), "tertiusDiff(net, nd$x)", "netterms(tertiusdiff(x))"),
             one("f", "f06", 6, "choice_coordination", list(kind = "und", n = 7, N = 50), "inertia(net) + indeg(net, ignoreRep = TRUE)", "inertia indeg ignorerep(indeg)"),
             one("f", "f07", 7, "choice_coordination", list(kind = "und", n = 7, N = 50), "inertia(net) + trans(net, ignoreRep = TRUE)", "inertia trans ignorerep(trans)"),
             one("f", "f08", 8, "choice_coordination", list(kind = "und", n = 8, N = 60), "inertia(net) + four(net, ignoreRep = TRUE)", "inertia four ignorerep(four)"),
             one("f", "f09", 9, "choice_coordination", list(kind = "und", n = 8, N = 60), "indeg(net) + nodeTrans(net, ignoreRep = TRUE) + tertius(net, nd$x, ignoreRep = TRUE)", "indeg nodetrans tertius(x) ignorerep(nodetrans tertius)"))
  sf2 <- randspecs("f", c(bf, pf), 40, 61, "choice_coordination", function(k) list(kind = "und", n = sample(6:9, 1), N = sample(40:70, 1), exo = EXOU, weights = (k %% 4 == 0),
         inttime = (k %% 5 == 0), ties = (k %% 3 != 0)), maxeff = 2, idoff = 10)
  specs <- c(specs, sf, sf2)
}
for (sp in specs) runspec(sp)
for (sec in secs) {
  r <- do.call(rbind, Filter(function(d) d$sec[1] == sec, res)); cd <- do.call(rbind, Filter(function(d) d$sec[1] == sec, cdf))
  if (!is.null(r)) { write.csv(r[, -1], file.path(out, paste0("results_", sec, ".csv")), row.names = FALSE); write.csv(cd[, -1], file.path(out, paste0("cases_", sec, ".csv")), row.names = FALSE) }
}
