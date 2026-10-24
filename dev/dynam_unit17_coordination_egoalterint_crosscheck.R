# DyNAM Part 0 / unit 17 crosscheck (the development notes) -
# choice_coordination effect expansion, "batch 14", resolved
# 2026-09-02 (continuing the goldfish-comparison list after batch 13's
# tertius()/four).
#
# Adds egoAlterInt() to choice_coordination's own effect set - the
# LAST remaining real-goldfish-eligible effect from batch 11's own
# eligibility sweep (dev/dynam_unit14_coordination_crosscheck.R's own
# header comment), completing that list (tie/tertiusDiff/mixedTrans
# remain a real, disclosed scope limit tied to the cross-network
# capability this engine doesn't have yet). Uses TWO distinct
# covariates (floor, floor2) - egoAlterInt on the SAME variable twice
# is not tested here since that degeneracy is already documented for
# plain choice and is not specific to coordination.

library(goldfish)

n <- 6
set.seed(42)
actors <- data.frame(label = paste0("A", 1:n), floor = c(1,1,2,2,3,3), floor2 = c(5,3,8,1,9,2))
nodes <- defineNodes(actors)
net <- defineNetwork(nodes = nodes, directed = FALSE)
sender <- sample(1:n, 12, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
events <- data.frame(time = 1:12, sender = paste0("A", sender),
                      receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

cat("=== choice_coordination ~ egoAlterInt(floor, floor2) ===\n")
m1 <- estimate(dep ~ egoAlterInt(nodes$floor, nodes$floor2), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m1)); print(logLik(m1))

cat("\n=== choice_coordination ~ inertia + egoAlterInt(floor, floor2) ===\n")
m2 <- estimate(dep ~ inertia(net) + egoAlterInt(nodes$floor, nodes$floor2), model = "DyNAM",
                subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m2)); print(logLik(m2))
