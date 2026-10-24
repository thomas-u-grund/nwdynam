# dev/dynam_unit21_rateintercept_twomode_crosscheck.R -- verifies the
# WITH-INTERCEPT rate sub-model combined with two-mode (bipartite)
# support, batch 18 (continuing the goldfish-parity list after batch 17's
# base intercept implementation). Real, unevenly-spaced timestamps
# (essential: the whole point of the intercept variant is sensitivity to
# actual elapsed time). Only `outdeg` is used as the covariate effect,
# since goldfish itself REJECTS `indeg` on a two-mode rate sub-model
# ("indeg effect must not use for type 'ego' when isTwoMode" - confirmed
# directly this session) - outdeg is the one rate-submodel effect already
# established as two-mode-eligible.
library(goldfish)
set.seed(7)
n1 <- 5  # mode 1 (senders)
n2 <- 4  # mode 2 (receivers)

actors1 <- data.frame(label = paste0("A", 1:n1), present = TRUE)
actors2 <- data.frame(label = paste0("B", 1:n2), present = TRUE)
nodes1 <- defineNodes(actors1)
nodes2 <- defineNodes(actors2)

mat0 <- matrix(0, n1, n2, dimnames = list(actors1$label, actors2$label))
twomodeNet <- defineNetwork(matrix = mat0, nodes = nodes1, nodes2 = nodes2, directed = TRUE)

# Every mode-1 actor must appear as sender at least once - otherwise
# `encode`-based reconstruction on the Stata side (which derives the
# actor universe purely from OBSERVED values, unlike goldfish's own
# explicit up-front actor declaration) would silently drop a real
# competing-risk actor from the hazard sum, producing a genuine mismatch
# that looks like a formula bug but isn't (found directly this way,
# 2026-09-02 - the very first version of this toy dataset happened not
# to include actor A1 as a sender in its own random 12-event sample).
nevents <- 12
repeat {
  events <- data.frame(
    time = cumsum(round(runif(nevents, 1, 15), 2)),
    sender = paste0("A", sample(1:n1, nevents, replace = TRUE)),
    receiver = paste0("B", sample(1:n2, nevents, replace = TRUE)),
    increment = 1
  )
  if (length(unique(events$sender)) == n1) break
}
twomodeNet <- linkEvents(x = twomodeNet, changeEvents = events, nodes = nodes1, nodes2 = nodes2)
depnet <- defineDependentEvents(events = events, nodes = nodes1, nodes2 = nodes2, defaultNetwork = twomodeNet)

mod <- estimate(depnet ~ 1 + outdeg, model = "DyNAM", subModel = "rate")
print(summary(mod))
cat(sprintf("INTERCEPT=%.8f OUTDEG=%.8f LOGLIK=%.8f\n",
    coef(mod)["Intercept"], coef(mod)["outdeg"], as.numeric(logLik(mod))))

write.csv(events, "dev/dynam_crosscheck_twomode_intercept_events.csv", row.names = FALSE)
cat(sprintf("N1=%d\n", n1))
