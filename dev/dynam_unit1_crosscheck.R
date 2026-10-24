# DyNAM Part 0 / unit 1 crosscheck (the development notes).
#
# Establishes the real reference implementation and a real, published
# crosscheck dataset - not simulated, unlike relevent's own bundled-nothing
# situation that forced nwrem's Part 0 into a simulate-then-recover design
# (the development notes). goldfish ships real datasets directly.
#
# Environment verified 2026-09-02: R 4.6.0 (aarch64-apple-darwin23),
# goldfish 1.6.12 installed from CRAN (install.packages("goldfish")).
# License: GPL (>= 3), confirmed directly via packageDescription("goldfish")
# - no unusual attribution requirements, unlike ergm's own
# GPL-3-plus-mandatory-statnet-attribution terms (the development notes).
#
# Dataset: Social_Evolution (bundled with goldfish) - the MIT "Social
# Evolution" college-dorm dataset, 84 actors, 439 timestamped phone-call
# events (`calls`), plus a separate `friendship` network and actor
# covariates (floor, gradeType). A real, previously-published dataset,
# not fabricated for this crosscheck.
#
# This script fits BOTH DyNAM sub-models separately (rate: who acts and
# when; choice: whom they choose) on identical data, capturing real
# reference coefficients for whatever from-scratch nwcommands
# implementation is eventually built against them. Re-run this script
# directly (not by re-typing the numbers below) when that implementation
# needs a live target to match - do not treat the numbers in this
# comment block as a substitute for actually running goldfish again.

library(goldfish)
data("Social_Evolution")

callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors,
                                         defaultNetwork = callNetwork)

cat("=== DyNAM choice sub-model: inertia + recip + indeg ===\n")
modChoice <- estimate(callsDependent ~ inertia + recip + indeg,
                       model = "DyNAM", subModel = "choice")
print(summary(modChoice))
# Reference (captured 2026-09-02, goldfish 1.6.12):
#   inertia  5.36556  (SE 0.19375)
#   recip    1.38944  (SE 0.16729)
#   indeg   -0.17353  (SE 0.09097)
#   logLik -695.46

cat("\n=== DyNAM rate sub-model: indeg + outdeg ===\n")
modRate <- estimate(callsDependent ~ indeg + outdeg,
                     model = "DyNAM", subModel = "rate")
print(summary(modRate))
# Reference (captured 2026-09-02, goldfish 1.6.12):
#   indeg   0.551445  (SE 0.066344)
#   outdeg  0.263784  (SE 0.028386)
#   logLik -1750.9
