# DyNAM Part 0 / unit 15 crosscheck (the development notes) -
# choice_coordination effect expansion, "batch 12", resolved
# 2026-09-02 (continuing the goldfish-comparison list after batch 11's
# choice_coordination core).
#
# Adds nodeTrans and trans to choice_coordination's own six v1
# effects. Both confirmed real-goldfish-eligible for
# choice_coordination in batch 11's own eligibility sweep
# (dev/dynam_unit14_coordination_crosscheck.R's own header comment).
# Mechanics turned out to need NO new full-matrix construction beyond
# what plain choice's own closure family already computes - see
# unw_dynam.do's own header comment (above dynam_coord_loglik_grad_
# multi()) for the exact reasoning, verified here on the SAME
# 6-actor/12-event toy network dev/dynam_unit14_coordination_
# crosscheck.R built (re-declared identically via the same seed).

library(goldfish)

n <- 6
set.seed(42)
actors <- data.frame(label = paste0("A", 1:n))
nodes <- defineNodes(actors)
net <- defineNetwork(nodes = nodes, directed = FALSE)
sender <- sample(1:n, 12, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
events <- data.frame(time = 1:12, sender = paste0("A", sender),
                      receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

cat("=== choice_coordination ~ nodeTrans ===\n")
m1 <- estimate(dep ~ nodeTrans(net), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m1)); print(logLik(m1))

cat("\n=== choice_coordination ~ trans ===\n")
m2 <- estimate(dep ~ trans(net), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m2)); print(logLik(m2))

cat("\n=== choice_coordination ~ inertia + nodeTrans + trans ===\n")
m3 <- estimate(dep ~ inertia(net) + nodeTrans(net) + trans(net), model = "DyNAM",
                subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m3)); print(logLik(m3))

write.csv(events, "dev/dynam_crosscheck_coord3_events.csv", row.names = FALSE)
