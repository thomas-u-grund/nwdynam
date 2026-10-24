# DyNAM Part 0 / unit 10 crosscheck (the development notes) - "effect
# expansion batch 7": tertius (both sub-models) and egoAlterInt
# (choice-only), continuing the systematic audit against goldfish's
# full effect catalog ("work your way through the list excluding
# dynami", 2026-09-02).
#
# tertius: default aggregateFun (mean) only, on the DEPENDENT network
# itself (goldfish's own default network argument) - v1 scope,
# matching every other structural effect in this package so far
# (custom aggregateFun/transformFun are not requestable).
#
# Real reference numbers below captured live against goldfish 1.6.12 on
# its own bundled Social_Evolution dataset (84 actors, 439 real
# phone-call events, own `floor`/`gradeType` covariates) - re-run this
# script directly, do not treat the numbers as a substitute for
# actually running goldfish again.

library(goldfish)
data("Social_Evolution")

callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors,
                                         defaultNetwork = callNetwork)

cat("=== choice ~ tertius(callNetwork, actors$floor) ===\n")
m1 <- estimate(callsDependent ~ tertius(callNetwork, actors$floor),
                model = "DyNAM", subModel = "choice")
print(summary(m1))

cat("\n=== choice ~ inertia + tertius(callNetwork, actors$floor) ===\n")
m2 <- estimate(callsDependent ~ inertia + tertius(callNetwork, actors$floor),
                model = "DyNAM", subModel = "choice")
print(summary(m2))

cat("\n=== rate ~ tertius(callNetwork, actors$floor) ===\n")
m3 <- estimate(callsDependent ~ tertius(callNetwork, actors$floor),
                model = "DyNAM", subModel = "rate")
print(summary(m3))

cat("\n=== choice ~ egoAlterInt(actors$floor, actors$gradeType) ===\n")
m4 <- estimate(callsDependent ~ egoAlterInt(actors$floor, actors$gradeType),
                model = "DyNAM", subModel = "choice")
print(summary(m4))

cat("\n=== choice ~ inertia + alter(actors$gradeType) + egoAlterInt(actors$floor, actors$gradeType) ===\n")
m5 <- estimate(callsDependent ~ inertia + alter(actors$gradeType) +
                  egoAlterInt(actors$floor, actors$gradeType),
                model = "DyNAM", subModel = "choice")
print(summary(m5))
