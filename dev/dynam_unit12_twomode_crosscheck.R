# DyNAM Part 0 / unit 12 crosscheck (the development notes) - two-mode
# (bipartite) DyNAM support, "expansion batch 9", resolved 2026-09-02,
# user: "continue" (after being shown the remaining large architectural
# items: two-mode, opportunitiesList, choice_coordination,
# cross-network effects, WITH-INTERCEPT rate).
#
# Toy two-mode affiliation network: 6 "people" (mode 1, node set
# `nodes`), 4 "orgs" (mode 2, node set `nodes2`), 60 events, people
# always sending to orgs, never the reverse - matching goldfish's own
# strictly one-directional two-mode DyNAM architecture (confirmed
# directly, not assumed - see unw_dynam.do's own header comment for
# the full investigation, including several effects the table
# advertises as `isTwoMode`-eligible that goldfish itself hard-rejects
# at runtime, and a deeper finding that even `same()`/`ego()` compute
# something different from a naive combined-covariate approach for
# two-mode, traced to a goldfish-internal row/column self-tie-index
# quirk not confidently reverse-engineered - so `nwdynam` ships ONLY
# `inertia`/`indeg` (choice)/`outdeg` (rate) for two-mode in this
# batch, verified below).
#
# Real reference numbers below captured live against goldfish 1.6.12 -
# re-run this script directly, do not treat the numbers as a
# substitute for actually running goldfish again.

library(goldfish)
set.seed(42)
n1 <- 6; n2 <- 4
shared <- c(3, 1, 4, 1, 5, 9, 2, 7, 1, 8)
people <- data.frame(label = paste0("P", 1:n1), sharedattr = shared[1:n1])
orgs <- data.frame(label = paste0("O", 1:n2), sharedattr = shared[(n1 + 1):(n1 + n2)])
nev <- 60
sender <- sample(1:n1, nev, replace = TRUE)
receiver <- sample(1:n2, nev, replace = TRUE)
events <- data.frame(time = 1:nev, sender = paste0("P", sender), receiver = paste0("O", receiver), increment = 1)
peopleNodes <- defineNodes(people)
orgNodes <- defineNodes(orgs)
affNet <- defineNetwork(nodes = peopleNodes, nodes2 = orgNodes, directed = TRUE)
affNet <- linkEvents(x = affNet, changeEvent = events, nodes = peopleNodes, nodes2 = orgNodes)
dep <- defineDependentEvents(events = events, nodes = peopleNodes, nodes2 = orgNodes, defaultNetwork = affNet)

cat("=== choice ~ inertia ===\n")
print(summary(estimate(dep ~ inertia, model = "DyNAM", subModel = "choice")))

cat("\n=== choice ~ indeg ===\n")
print(summary(estimate(dep ~ indeg, model = "DyNAM", subModel = "choice")))

cat("\n=== choice ~ inertia + indeg ===\n")
print(summary(estimate(dep ~ inertia + indeg, model = "DyNAM", subModel = "choice")))

cat("\n=== rate ~ outdeg ===\n")
print(summary(estimate(dep ~ outdeg, model = "DyNAM", subModel = "rate")))

# Write out the exact event data used, for the direct-Mata-call .do
# crosscheck and the cscripts/test_nwdynam_twomode_ado.do .ado-level
# test to both build on the SAME real dataset.
write.csv(events, "dev/dynam_crosscheck_twomode_events.csv", row.names = FALSE)
write.csv(data.frame(id = 1:(n1 + n2), shared = shared), "dev/dynam_crosscheck_twomode_actorcov.csv", row.names = FALSE)
cat("\nn1=", n1, " n2=", n2, "\n")
