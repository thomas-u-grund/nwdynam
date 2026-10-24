# DyNAM Part 0 / unit 13 crosscheck (the development notes) -
# opportunitiesList support, "expansion batch 10", resolved 2026-09-02,
# user: "work your way through it" (continuing the goldfish-comparison
# list after two-mode support).
#
# goldfish's own `estimate()` argument `opportunitiesList` (choice
# sub-model only, per its own documentation: "ONLY for choice models")
# restricts the candidate risk set PER EVENT to a caller-supplied list
# of "available" nodes. Mechanics derived from a real toy example
# FIRST, not guessed: `opportunitiesList` takes NUMERIC node indices
# (1-based, matching each actor's own row position in the `nodes` data
# frame) - passing character labels instead silently produces a
# confusing "Active node k not available in event i" error, even when
# the label IS logically present, because goldfish does not coerce
# labels to indices for this argument. The sender's own self-exclusion
# still applies unconditionally regardless of whether the sender is
# itself listed in the per-event opportunity set. Confirmed as a
# genuine restriction (not a no-op) by comparing a full-risk-set
# opportunitiesList (matches the no-opportunitiesList baseline exactly)
# against one that excludes an actor from every event (changes both
# the coefficient and the log-likelihood).
#
# Real reference numbers below captured live against goldfish 1.6.12 -
# re-run this script directly, do not treat the numbers as a
# substitute for actually running goldfish again.

library(goldfish)
set.seed(7)
n <- 8
actors <- data.frame(label = paste0("A", 1:n), attr1 = c(3, 1, 4, 1, 5, 9, 2, 6))
nodes <- defineNodes(actors)
net <- defineNetwork(nodes = nodes, directed = TRUE)
nev <- 30
sender <- sample(1:n, nev, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
events <- data.frame(time = 1:nev, sender = paste0("A", sender), receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

cat("=== baseline: choice ~ inertia + indeg, no opportunitiesList ===\n")
m0 <- estimate(dep ~ inertia + indeg, model = "DyNAM", subModel = "choice")
print(summary(m0))

# Genuinely different opportunity sets per event: exclude ONE random
# actor (never the sender or the true receiver) from each event's own
# risk set.
set.seed(99)
oppList <- lapply(1:nev, function(i) {
  s <- sender[i]; r <- receiver[i]
  candidates <- setdiff(1:n, s)
  excl_pool <- setdiff(candidates, r)
  excl <- sample(excl_pool, 1)
  setdiff(1:n, excl)
})

cat("\n=== choice ~ inertia + indeg, WITH per-event opportunitiesList (1 random exclusion each) ===\n")
m1 <- estimate(dep ~ inertia + indeg, model = "DyNAM", subModel = "choice",
                estimationInit = list(opportunitiesList = oppList))
print(summary(m1))

# Regression check: full risk set every event reproduces the baseline exactly.
oppListFull <- lapply(1:nev, function(i) 1:n)
cat("\n=== choice ~ inertia + indeg, opportunitiesList = full risk set (must match baseline exactly) ===\n")
m2 <- estimate(dep ~ inertia + indeg, model = "DyNAM", subModel = "choice",
                estimationInit = list(opportunitiesList = oppListFull))
print(summary(m2))

write.csv(events, "dev/dynam_crosscheck_oppl_events.csv", row.names = FALSE)
oppdf <- do.call(rbind, lapply(seq_along(oppList), function(i) data.frame(event = i, actor = oppList[[i]])))
write.csv(oppdf, "dev/dynam_crosscheck_oppl_availability.csv", row.names = FALSE)
