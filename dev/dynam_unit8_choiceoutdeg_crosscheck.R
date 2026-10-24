# DyNAM Part 0 / unit 8 crosscheck (the development notes) - adds the
# choice sub-model's own `outdeg` effect (goldfish's "alter type"
# outdegree), correcting nwdynam.sthlp's own "outdeg: rate sub-model
# only" label, which was checked against goldfish's own
# doc/goldfishEffects.Rmd and found WRONG: the "Structural Effects"
# table there lists outdeg as valid for BOTH rate and choice, not
# rate-only - a genuine missing effect, not a correct scope
# restriction (user report, 2026-09-02).
#
# Real reference numbers below captured live against goldfish 1.6.12 on
# its own bundled Social_Evolution dataset (84 actors, 439 real
# phone-call events) - re-run this script directly, do not treat the
# numbers as a substitute for actually running goldfish again.

library(goldfish)
data("Social_Evolution")

callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors,
                                         defaultNetwork = callNetwork)

cat("=== choice ~ outdeg (alter type - candidate receiver's own out-degree) ===\n")
m0 <- estimate(callsDependent ~ outdeg, model = "DyNAM", subModel = "choice")
print(summary(m0))

cat("\n=== choice ~ inertia + outdeg ===\n")
m1 <- estimate(callsDependent ~ inertia + outdeg, model = "DyNAM", subModel = "choice")
print(summary(m1))

cat("\n=== choice ~ inertia + recip + indeg + outdeg (all four structural effects together) ===\n")
m2 <- estimate(callsDependent ~ inertia + recip + indeg + outdeg, model = "DyNAM", subModel = "choice")
print(summary(m2))
