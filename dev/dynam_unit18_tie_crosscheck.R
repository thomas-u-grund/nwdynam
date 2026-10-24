# DyNAM Part 0 / unit 18 crosscheck (the development notes) -
# tie() support ("cross-network effects, v1 scope"), resolved
# 2026-09-02, continuing the goldfish-comparison list after batch 14's
# choice_coordination egoAlterInt().
#
# goldfish's own `tie(network, ...)` effect: `s(i,j,t,x) = I(x_ij>0)`
# (unweighted) - identical in shape to `inertia()` but reading a
# SEPARATE, exogenous network `x` instead of the dependent network.
# v1 scope (a real, disclosed limit, not silently narrow): the
# exogenous network here is STATIC (declared via `defineNetwork(matrix
# = ..., nodes = nodes)`, no `changeEvent` linked, so it never changes
# across events) - genuinely dynamically-evolving cross-network
# effects (where the SECOND network is itself event-declared) remain a
# real, disclosed follow-on for a later batch.

library(goldfish)

n <- 6
set.seed(42)
actors <- data.frame(label = paste0("A", 1:n))
nodes <- defineNodes(actors)

net <- defineNetwork(nodes = nodes, directed = TRUE)
sender <- sample(1:n, 12, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
events <- data.frame(time = 1:12, sender = paste0("A", sender),
                      receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

set.seed(7)
staticmat <- matrix(sample(0:1, n * n, replace = TRUE, prob = c(0.6, 0.4)), n, n)
diag(staticmat) <- 0
staticnet <- defineNetwork(matrix = staticmat, nodes = nodes, directed = TRUE)

cat("=== choice ~ tie(staticnet) ===\n")
m1 <- estimate(dep ~ tie(staticnet), model = "DyNAM", subModel = "choice",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m1)); print(logLik(m1))

cat("\n=== choice ~ inertia + tie(staticnet) ===\n")
m2 <- estimate(dep ~ inertia(net) + tie(staticnet), model = "DyNAM", subModel = "choice",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m2)); print(logLik(m2))

write.csv(events, "dev/dynam_crosscheck_tie_events.csv", row.names = FALSE)
write.csv(as.data.frame(staticmat), "dev/dynam_crosscheck_tie_staticmat.csv", row.names = FALSE)
