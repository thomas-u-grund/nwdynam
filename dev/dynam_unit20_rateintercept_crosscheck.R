# DyNAM Part 0 / unit 20 crosscheck (the development notes) -
# WITH-INTERCEPT rate sub-model, "expansion batch 17", resolved
# 2026-09-02, continuing the goldfish-comparison list after
# tertiusDiff's own investigation (found genuinely too deep to
# reverse-engineer safely, disclosed rather than guessed).
#
# The genuinely continuous-time competing-risks hazard variant this
# project's own docs always described as deferred - unlike the
# NO-intercept rate model (which reduces to the same ordinal partial
# likelihood as everything else in this package), an intercept makes
# the model sensitive to the ACTUAL SCALE of elapsed real time between
# events. Requested via a literal `1` as the first formula term
# (confirmed directly via goldfish's own `parseIntercept()` source);
# `dep ~ 1` ALONE (no other term) hits a real goldfish-internal bug
# unrelated to this package (`strsplit(): non-character argument`) -
# not investigated further, worked around by always including at
# least one real effect alongside the intercept.

library(goldfish)

n <- 6
set.seed(42)
actors <- data.frame(label = paste0("A", 1:n))
nodes <- defineNodes(actors)
net <- defineNetwork(nodes = nodes, directed = TRUE)
sender <- sample(1:n, 10, replace = TRUE)
receiver <- sapply(sender, function(s) sample(setdiff(1:n, s), 1))
# REAL, unevenly-spaced timestamps (not sequential integers) - the
# genuinely time-sensitive variant needs real elapsed time to matter.
times <- cumsum(round(runif(10, 1, 20), 2))
events <- data.frame(time = times, sender = paste0("A", sender),
                      receiver = paste0("A", receiver), increment = 1)
net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

cat("=== rate ~ 1 + indeg(net) ===\n")
m1 <- estimate(dep ~ 1 + indeg(net), model = "DyNAM", subModel = "rate")
print(coef(m1)); print(logLik(m1))

cat("\n=== rate ~ 1 + outdeg(net) ===\n")
m2 <- estimate(dep ~ 1 + outdeg(net), model = "DyNAM", subModel = "rate")
print(coef(m2)); print(logLik(m2))

write.csv(events, "dev/dynam_crosscheck_intercept_events.csv", row.names = FALSE)
