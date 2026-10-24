# DyNAM Part 0 / unit 7 crosscheck (the development notes) - extends
# window() to the RATE sub-model's indeg/outdeg effects, per the user's
# follow-up request "can we make [window] also available for the rate
# part of the model?" (2026-09-02, immediately after the choice-side
# inertiawindow()/recipwindow() self-activation redesign).
#
# goldfish's own doc/goldfishEffects.Rmd effect table (row "indeg"/
# "outdeg" under "Node/actor statistics") lists `window` as a real
# argument on BOTH structural node effects, and its own
# "Structural Effects" table confirms indeg/outdeg apply to BOTH the
# rate and choice sub-models (indeg: rate=yes, choice=yes; outdeg:
# rate=yes, choice=no under choice_coordination but yes under plain
# choice/rate). This crosscheck targets the RATE sub-model only, since
# that is the scope gap the development notes' own "What is explicitly
# NOT yet implemented" section already flagged ("a windowed rate-
# sub-model effect is also not yet wired").
#
# Real reference numbers below captured live against goldfish 1.6.12 on
# its own bundled Social_Evolution dataset (84 actors, 439 real phone-call
# events) - re-run this script directly, do not treat the numbers as a
# substitute for actually running goldfish again.

library(goldfish)
data("Social_Evolution")

win7d <- 7 * 24 * 60 * 60
win1d <- 1 * 24 * 60 * 60

callNetwork <- defineNetwork(nodes = actors, directed = TRUE)
callNetwork <- linkEvents(x = callNetwork, changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(events = calls, nodes = actors,
                                         defaultNetwork = callNetwork)

cat("=== rate ~ indeg, no window (regression check against unit1's own numbers) ===\n")
m0 <- estimate(callsDependent ~ indeg, model = "DyNAM", subModel = "rate")
print(summary(m0))

cat("\n=== rate ~ indeg, window = 7 days ===\n")
m1 <- estimate(callsDependent ~ indeg(callNetwork, window = win7d),
                model = "DyNAM", subModel = "rate")
print(summary(m1))

cat("\n=== rate ~ indeg, window = 1 day ===\n")
m2 <- estimate(callsDependent ~ indeg(callNetwork, window = win1d),
                model = "DyNAM", subModel = "rate")
print(summary(m2))

cat("\n=== rate ~ outdeg, window = 1 day ===\n")
m3 <- estimate(callsDependent ~ outdeg(callNetwork, window = win1d),
                model = "DyNAM", subModel = "rate")
print(summary(m3))

cat("\n=== rate ~ indeg(window=7 days) + outdeg(window=1 day) - GENUINELY DIFFERENT windows on both simultaneously ===\n")
m4 <- estimate(callsDependent ~ indeg(callNetwork, window = win7d) +
                  outdeg(callNetwork, window = win1d),
                model = "DyNAM", subModel = "rate")
print(summary(m4))
