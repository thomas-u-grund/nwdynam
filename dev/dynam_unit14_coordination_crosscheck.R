# DyNAM Part 0 / unit 14 crosscheck (the development notes) -
# choice_coordination support, "expansion batch 11", resolved
# 2026-09-02 (user: "do choice_coordination now").
#
# goldfish's own third DyNAM sub-model: a "multinomial-multinomial"
# joint model over UNORDERED actor pairs, for undirected
# ("coordination") tie-formation events - a genuinely different
# likelihood from `choice` (Stadtfeld, Hollway & Block 2017, "Dynamic
# Network Actor Models: Investigating Coordination Ties through Time,"
# Sociological Methodology 47(1), equations 6/9/11). The formula was
# derived by reading goldfish's own C++ source directly
# (src/compute_coordination_selection.cpp, from the CRAN 1.6.12 source
# tarball - NOT the paper's own prose, and NOT the vignette's own
# effect-eligibility table, which self-contradicts on `alter`), then
# independently reconstructed and verified in R against the real
# goldfish fits below to 5-6 decimal places - see unw_dynam.do's own
# header comment (above dynam_coord_loglik_grad_multi()) for the full
# formula and its gradient.
#
# Effect eligibility for choice_coordination was checked DIRECTLY
# against real goldfish (preprocessingOnly=TRUE probes for every
# candidate effect name), not trusted from goldfish's own
# goldfishEffects.Rmd vignette table (already known unreliable in both
# directions from the two-mode batch, and directly caught contradicting
# its own prose text on `alter` here too).

library(goldfish)

## --- Toy 1: 4 actors, 4 events, inertia only ---
n <- 4
actors1 <- data.frame(label = paste0("A", 1:n))
nodes1 <- defineNodes(actors1)
net1 <- defineNetwork(nodes = nodes1, directed = FALSE)
events1 <- data.frame(
  time = 1:4,
  sender = c("A1", "A2", "A1", "A3"),
  receiver = c("A2", "A3", "A3", "A4"),
  increment = 1
)
net1 <- linkEvents(x = net1, changeEvent = events1, nodes = nodes1)
dep1 <- defineDependentEvents(events = events1, nodes = nodes1, defaultNetwork = net1)

cat("=== Toy 1: choice_coordination ~ inertia (4 actors, 4 events) ===\n")
m1 <- estimate(dep1 ~ inertia, model = "DyNAM", subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m1))
cat("coef:\n"); print(coef(m1))
cat("logLik:\n"); print(logLik(m1))

write.csv(events1, "dev/dynam_crosscheck_coord1_events.csv", row.names = FALSE)

## --- Toy 2: 6 actors, 12 events, inertia + same(floor) ---
n2 <- 6
set.seed(42)
actors2 <- data.frame(label = paste0("A", 1:n2), floor = c(1,1,2,2,3,3))
nodes2 <- defineNodes(actors2)
net2 <- defineNetwork(nodes = nodes2, directed = FALSE)
sender2 <- sample(1:n2, 12, replace = TRUE)
receiver2 <- sapply(sender2, function(s) sample(setdiff(1:n2, s), 1))
events2 <- data.frame(time = 1:12, sender = paste0("A", sender2),
                       receiver = paste0("A", receiver2), increment = 1)
net2 <- linkEvents(x = net2, changeEvent = events2, nodes = nodes2)
dep2 <- defineDependentEvents(events = events2, nodes = nodes2, defaultNetwork = net2)

cat("\n=== Toy 2: choice_coordination ~ inertia + same(floor) (6 actors, 12 events) ===\n")
m2 <- estimate(dep2 ~ inertia + same(nodes2$floor), model = "DyNAM",
                subModel = "choice_coordination",
                estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m2))
cat("coef:\n"); print(coef(m2))
cat("logLik:\n"); print(logLik(m2))

write.csv(events2, "dev/dynam_crosscheck_coord2_events.csv", row.names = FALSE)
write.csv(actors2, "dev/dynam_crosscheck_coord2_actors.csv", row.names = FALSE)

## --- Toy 3: same 6-actor data, the remaining v1 effects (indeg,
## diff, sim, alter) - each FIT ALONE, not combined: diff()/sim() on
## the SAME covariate are exactly negatives of each other (same
## degeneracy already documented for plain choice), and combining all
## four on this small a toy network hit a genuine "matrix cannot be
## inverted" collinearity failure when first drafted - not a formula
## bug, the identical finding this session already made for plain
## choice's diff()/sim() (see unw_dynam.do's own "expansion batch 1"
## comment) applying here too. Fit separately instead, matching
## dev/dynam_unit4_attreffects_crosscheck.R's own precedent for the
## same reason.
cat("\n=== Toy 3a: choice_coordination ~ indeg ===\n")
m3a <- estimate(dep2 ~ indeg(net2), model = "DyNAM", subModel = "choice_coordination",
                 estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m3a))
cat("coef:\n"); print(coef(m3a))
cat("logLik:\n"); print(logLik(m3a))

cat("\n=== Toy 3b: choice_coordination ~ diff(floor) ===\n")
m3b <- estimate(dep2 ~ diff(nodes2$floor), model = "DyNAM", subModel = "choice_coordination",
                 estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m3b))
cat("coef:\n"); print(coef(m3b))
cat("logLik:\n"); print(logLik(m3b))

cat("\n=== Toy 3c: choice_coordination ~ sim(floor) ===\n")
m3c <- estimate(dep2 ~ sim(nodes2$floor), model = "DyNAM", subModel = "choice_coordination",
                 estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m3c))
cat("coef:\n"); print(coef(m3c))
cat("logLik:\n"); print(logLik(m3c))

cat("\n=== Toy 3d: choice_coordination ~ alter(floor) ===\n")
m3d <- estimate(dep2 ~ alter(nodes2$floor), model = "DyNAM", subModel = "choice_coordination",
                 estimationInit = list(maxIterations = 40, initialDamping = 10))
print(summary(m3d))
cat("coef:\n"); print(coef(m3d))
cat("logLik:\n"); print(logLik(m3d))
