# DyNAM Unit 2 (rate sub-model) formula derivation - the key empirical
# finding that CORRECTED unw_dynam.do's own original Unit 1 header
# comment, which had characterized goldfish's ENTIRE rate sub-model as a
# genuine continuous-time competing-risks hazard on real inter-event
# waiting times. That is only true of the WITH-INTERCEPT variant
# (goldfish's own teaching1.Rmd vignette: "an intercept of -14 means a
# waiting time of 334 hours" - a real-time-scale-sensitive statement).
#
# This script fits a tiny, fully hand-traceable 3-actor/4-event toy
# example with goldfish's own rate sub-model, WITHOUT an intercept, and
# records the exact logLik - used to numerically reverse-engineer the
# underlying formula (see dev/dynam_unit2_rate_crosscheck.do for the
# from-scratch Mata engine, dynam_rate_loglik_grad_unit1() in
# unw_dynam.do, verified against BOTH this toy example and the real
# Social_Evolution dataset).
#
# Result (captured 2026-09-02): the no-intercept rate sub-model's
# logLik matches a pure ORDINAL PARTIAL LIKELIHOOD (Cox-style, no real
# time values used at all - only event ORDER) over the full n-actor
# sender risk set, using ego-type indeg/outdeg (an actor's own
# in/out-degree - count of distinct partners, unweighted) as predictors -
# to full double precision, both on this toy example and on the real
# Social_Evolution dataset (dev/dynam_unit1_crosscheck.R's own reference
# numbers: indeg=0.551445, outdeg=0.263784, logLik=-1750.9 - re-derived
# here from a from-scratch R implementation of the hypothesized ordinal
# formula, matching to 3+ decimal places and recovering identical MLE
# coefficients from a fresh optimize() run).

library(goldfish)

actors <- data.frame(label = c("A1","A2","A3"), present = TRUE)
calls <- data.frame(
  time = c(1, 3, 8, 20),
  sender = c("A1","A2","A1","A3"),
  receiver = c("A2","A1","A3","A1"),
  increment = 1
)

net <- defineNetwork(nodes = actors, directed = TRUE)
net <- linkEvents(x = net, changeEvent = calls, nodes = actors)
dep <- defineDependentEvents(events = calls, nodes = actors, defaultNetwork = net)

cat("=== goldfish rate submodel, NO intercept, outdeg(net) only ===\n")
m1 <- estimate(dep ~ outdeg(net), model = "DyNAM", subModel = "rate")
print(summary(m1))
# Reference (captured 2026-09-02): outdeg=-1.0050, logLik=-3.9196

cat("\n=== goldfish rate submodel, WITH intercept, outdeg(net) only (NOT matched by the ordinal ===\n")
cat("=== formula below - this is the genuinely-different, still-deferred continuous-time variant) ===\n")
m2 <- estimate(dep ~ 1 + outdeg(net), model = "DyNAM", subModel = "rate")
print(summary(m2))
# Reference (captured 2026-09-02): Intercept=-1.91522, outdeg=-1.46981, logLik=-13.131

# From-scratch R re-implementation of the hypothesized ordinal formula -
# independent of goldfish's own estimate() internals, used to confirm
# the match is not a coincidence of this one toy example.
n <- 3
sid <- match(calls$sender, actors$label)
rid <- match(calls$receiver, actors$label)
nevents <- nrow(calls)

loglik_ordinal_rate <- function(theta) {
  tie <- matrix(0, n, n)
  ll <- 0
  for (k in 1:nevents) {
    s <- sid[k]; r <- rid[k]
    outdeg_ego <- rowSums(tie > 0)
    lp <- theta[1] * outdeg_ego
    mx <- max(lp)
    denom <- mx + log(sum(exp(lp - mx)))
    ll <- ll + (lp[s] - denom)
    tie[s, r] <- 1
  }
  ll
}

cat("\n=== independent from-scratch ordinal-partial-likelihood re-derivation ===\n")
cat("logLik at goldfish's own reference theta=-1.0050:", loglik_ordinal_rate(-1.0050), "\n")
cat("(goldfish's own no-intercept reference logLik: -3.9196)\n")

opt <- optim(0, function(th) -loglik_ordinal_rate(th), method = "BFGS", hessian = TRUE)
cat("recovered theta_hat:", opt$par, "  logLik:", -opt$value, "\n")
cat("(matches goldfish's own no-intercept fit: theta=-1.0050, logLik=-3.9196)\n")
