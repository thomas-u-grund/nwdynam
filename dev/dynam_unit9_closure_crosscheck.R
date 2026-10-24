# DyNAM Part 0 / unit 9 crosscheck (the development notes) - "effect
# expansion batch 6": nodeTrans (both sub-models) and the choice-only
# two-path closure effects trans/cycle/commonSender/commonReceiver/four,
# all evaluated on the DEPENDENT network itself (goldfish's own default
# when no explicit network argument is given - "For network effects,
# this would be either an explanatory network or the dependent network
# (the default)," doc/goldfishEffects.Rmd's own "Effect arguments"
# section). Exogenous/cross-network use of these same effects (an
# explicit different network argument), and the two-network
# mixedTrans/mixedCycle/mixedCommonSender/mixedCommonReceiver variants,
# are a separate, larger architectural extension (nwdynam has no
# mechanism yet for an effect to read a network OTHER than the one
# being predicted) and are intentionally out of scope for this batch -
# disclosed directly, not silently narrowed.
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

cat("=== choice ~ trans ===\n")
m1 <- estimate(callsDependent ~ trans, model = "DyNAM", subModel = "choice")
print(summary(m1))

cat("\n=== choice ~ inertia + trans ===\n")
m2 <- estimate(callsDependent ~ inertia + trans, model = "DyNAM", subModel = "choice")
print(summary(m2))

cat("\n=== choice ~ cycle ===\n")
m3 <- estimate(callsDependent ~ cycle, model = "DyNAM", subModel = "choice")
print(summary(m3))

cat("\n=== choice ~ commonSender ===\n")
m4 <- estimate(callsDependent ~ commonSender, model = "DyNAM", subModel = "choice")
print(summary(m4))

cat("\n=== choice ~ commonReceiver ===\n")
m5 <- estimate(callsDependent ~ commonReceiver, model = "DyNAM", subModel = "choice")
print(summary(m5))

cat("\n=== choice ~ four ===\n")
m6 <- estimate(callsDependent ~ four, model = "DyNAM", subModel = "choice")
print(summary(m6))

cat("\n=== choice ~ nodeTrans ===\n")
m7 <- estimate(callsDependent ~ nodeTrans, model = "DyNAM", subModel = "choice")
print(summary(m7))

cat("\n=== rate ~ nodeTrans ===\n")
m8 <- estimate(callsDependent ~ nodeTrans, model = "DyNAM", subModel = "rate")
print(summary(m8))

cat("\n=== choice ~ inertia + recip + trans + cycle + commonSender + commonReceiver (all together) ===\n")
m9 <- estimate(callsDependent ~ inertia + recip + trans + cycle + commonSender + commonReceiver,
                model = "DyNAM", subModel = "choice")
print(summary(m9))
