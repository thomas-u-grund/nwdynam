# nwdynam: events with identical times (ties) -- goldfish 1.6.12 reference numbers for cscripts/test_nwdynam_ties.do.
# goldfish treats events with the same time as SIMULTANEOUS: the statistics of every event of a time point are computed
# from the state before that time point (checked: the fit does not change when the tied rows are permuted).
# Usage: Rscript dev/dynam_ties_crosscheck.R <path to socialevolution.dta of the core package>
# Writes dev/dynam_ties_toy_events.csv (the toy data used in the test) and prints all numbers.
suppressMessages({library(goldfish); library(haven)})
args <- commandArgs(TRUE)
d <- read_dta(args[1])
actors <- data.frame(label = sort(unique(c(d$sender, d$receiver))), present = TRUE, stringsAsFactors = FALSE)
nodes <- defineNodes(actors)
fit <- function(ev, form, sub, nodes, directed = TRUE) {
  net <- defineNetwork(nodes = nodes, directed = directed)
  net <- linkEvents(x = net, changeEvent = ev, nodes = nodes)
  dep <- defineDependentEvents(events = ev, nodes = nodes, defaultNetwork = net)
  assign("dep", dep, envir = globalenv()); assign("net", net, envir = globalenv()); assign("ev", ev, envir = globalenv()); assign("nodes", nodes, envir = globalenv())
  m <- suppressWarnings(estimate(as.formula(paste("dep ~", form)), model = "DyNAM", subModel = sub,
                                 estimationInit = list(maxIterations = 60, initialDamping = 10)))
  r <- c(coef(m), ll = as.numeric(m$logLikelihood))
  cat(sprintf("%-14s %s\n", sub, form)); print(signif(r, 9)); invisible(r)
}
mk <- function(f) { ev <- data.frame(time = f(d$t), sender = d$sender, receiver = d$receiver, increment = 1,
                                     stringsAsFactors = FALSE); ev[order(ev$time), ] }
cat("##### A: socialevolution, time = ceil(t/20)  (4 time points)\n")
evA <- mk(function(t) ceiling(t / 20))
fit(evA, "inertia + recip + indeg", "choice", nodes)
fit(evA, "indeg + outdeg", "rate", nodes)
cat("##### B: socialevolution, time = ceil(2t)  (47 time points)\n")
evB <- mk(function(t) ceiling(t * 2))
fit(evB, "inertia + recip + indeg", "choice", nodes)
fit(evB, "inertia + recip + indeg + outdeg(net) + trans(net)", "choice", nodes)
fit(evB, "inertia(net, window = 5) + recip + indeg", "choice", nodes)
fit(evB, "indeg + outdeg", "rate", nodes)
fit(evB, "1 + indeg(net) + outdeg(net)", "rate", nodes)
# permuting the rows of tied events does not change goldfish's result
set.seed(7); evB2 <- evB[sample(nrow(evB)), ]; evB2 <- evB2[order(evB2$time), ]
cat("##### B with the tied rows permuted\n")
fit(evB2, "inertia + recip + indeg", "choice", nodes)
cat("##### C: toy undirected data with tied times, choice_coordination\n")
set.seed(3)
n <- 6
toyact <- data.frame(label = paste0("A", 1:n))
tnodes <- defineNodes(toyact)
s <- sample(1:n, 14, replace = TRUE)
r <- sapply(s, function(x) sample(setdiff(1:n, x), 1))
evC <- data.frame(time = rep(1:7, each = 2), sender = paste0("A", s), receiver = paste0("A", r), increment = 1,
                  stringsAsFactors = FALSE)
write.csv(evC, "dev/dynam_ties_toy_events.csv", row.names = FALSE)
fit(evC, "inertia", "choice_coordination", tnodes, directed = FALSE)
cat("##### D: the same toy events, directed\n")
fit(evC, "inertia + recip + indeg", "choice", tnodes)
fit(evC, "indeg + outdeg", "rate", tnodes)
fit(evC, "1 + indeg(net) + outdeg(net)", "rate", tnodes)
