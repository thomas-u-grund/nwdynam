library(goldfish)
data("Social_Evolution")
callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors, defaultNetwork = callNetwork)

cat("=== choice ~ inertia(weighted=TRUE) ===\n")
m1 <- estimate(callsDependent ~ inertia(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "choice")
print(summary(m1))

cat("\n=== choice ~ recip(weighted=TRUE) ===\n")
m2 <- estimate(callsDependent ~ recip(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "choice")
print(summary(m2))

cat("\n=== choice ~ indeg(weighted=TRUE) ===\n")
m3 <- estimate(callsDependent ~ indeg(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "choice")
print(summary(m3))

cat("\n=== rate ~ outdeg(weighted=TRUE) ===\n")
m4 <- estimate(callsDependent ~ outdeg(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "rate")
print(summary(m4))

cat("\n=== choice ~ inertia(weighted=TRUE) + recip(weighted=TRUE) ===\n")
m5 <- estimate(callsDependent ~ inertia(callNetwork, weighted = TRUE) + recip(callNetwork, weighted = TRUE),
                model = "DyNAM", subModel = "choice")
print(summary(m5))
library(goldfish)
data("Social_Evolution")
callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors, defaultNetwork = callNetwork)

cat("=== choice ~ outdeg(weighted=TRUE) ===\n")
m1 <- estimate(callsDependent ~ outdeg(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "choice")
print(summary(m1))

cat("\n=== rate ~ indeg(weighted=TRUE) ===\n")
m2 <- estimate(callsDependent ~ indeg(callNetwork, weighted = TRUE), model = "DyNAM", subModel = "rate")
print(summary(m2))
