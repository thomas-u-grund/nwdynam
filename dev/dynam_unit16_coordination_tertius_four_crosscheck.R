# DyNAM Part 0 / unit 16 crosscheck (the development notes) -
# choice_coordination effect expansion, "batch 13", resolved
# 2026-09-02 (continuing the goldfish-comparison list after batch 12's
# nodeTrans/trans).
#
# Adds tertius() and four to choice_coordination's own effect set.
# Both confirmed real-goldfish-eligible for choice_coordination in
# batch 11's own eligibility sweep (dev/dynam_unit14_coordination_
# crosscheck.R's own header comment). `four` needed a genuinely NEW
# full-matrix adaptation of plain choice's own distinct-actors
# correction, with the diagonal correction term indexed by the
# CANDIDATE, not the sender - a real transcription pitfall caught
# directly by re-deriving in R first (a first draft indexed it by the
# sender instead and reproduced a real but WRONG number, -32.039 vs.
# goldfish's own -32.288). `tertius` needed no new construction at all
# - see unw_dynam.do's own header comment for the exact reasoning,
# verified here on the SAME 6-actor/12-event toy network dev/
# dynam_unit14_coordination_crosscheck.R built (re-declared identically
# via the same seed).

library(goldfish)

n <- 6
set.seed(42)
actors <- data.frame(label = paste0("A", 1:n), floor = c(1,1,2,2,3,3))
nodes <- defineNodes(actors)
net <- defineNetwork(nodes = nodes, directed = FALSE)
sender <- sample(1:n, 12, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
events <- data.frame(time = 1:12, sender = paste0("A", sender),
                      receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

cat("=== choice_coordination ~ tertius(floor) ===\n")
m1 <- estimate(dep ~ tertius(net, nodes$floor), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m1)); print(logLik(m1))

cat("\n=== choice_coordination ~ four ===\n")
m2 <- estimate(dep ~ four(net), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m2)); print(logLik(m2))

cat("\n=== choice_coordination ~ inertia + tertius(floor) + four ===\n")
m3 <- estimate(dep ~ inertia(net) + tertius(net, nodes$floor) + four(net), model = "DyNAM",
                subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m3)); print(logLik(m3))
