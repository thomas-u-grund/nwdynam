# DyNAM Part 0 / unit 19 crosscheck (the development notes) -
# choice_coordination effect expansion, "batch 16", resolved
# 2026-09-02, continuing the goldfish-comparison list after batch 15's
# tie() for submodel(choice).
#
# Extends tie() (cross-network effects, v1 scope: a SEPARATE, STATIC
# exogenous network) to choice_coordination. Mechanically the simplest
# possible addition: unlike plain choice's own `tie_i = tiemat[s,.]'
# (one row), choice_coordination already needs the FULL n x n
# statistic matrix for every effect, so `S_tie' is simply `tiemat'
# itself - no computation at all. Verified here on a fresh 6-actor/
# 12-event undirected toy network with a genuinely separate static
# exogenous network, both alone and combined with `inertia`.

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

set.seed(7)
staticmat <- matrix(sample(0:1, n * n, replace = TRUE, prob = c(0.6, 0.4)), n, n)
diag(staticmat) <- 0
staticnet <- defineNetwork(matrix = staticmat, nodes = nodes, directed = TRUE)

cat("=== choice_coordination ~ tie(staticnet) ===\n")
m1 <- estimate(dep ~ tie(staticnet), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m1)); print(logLik(m1))

cat("\n=== choice_coordination ~ inertia + tie(staticnet) ===\n")
m2 <- estimate(dep ~ inertia(net) + tie(staticnet), model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(coef(m2)); print(logLik(m2))

write.csv(events, "dev/dynam_crosscheck_coordtie_events.csv", row.names = FALSE)
write.csv(as.data.frame(staticmat), "dev/dynam_crosscheck_coordtie_staticmat.csv", row.names = FALSE)
