# nwdynam 1.8.1: goldfish 1.6.12 reference for test_nwdynam_fix181.do (first 631 emails, reporting hierarchy as an event network at time 0 and as a static network).
# Usage: Rscript dev/dynam_fix181_crosscheck.R  (run from the repository root)
suppressMessages(library(goldfish))
ea <- read.csv("dev/fix181_actors.csv")$label
em <- read.csv("dev/fix181_emails.csv"); hi <- read.csv("dev/fix181_hier.csv")
nodes <- defineNodes(data.frame(label = ea))
mail <- defineNetwork(nodes = nodes, directed = TRUE)
events <- data.frame(time = em$t, sender = as.character(em$sender), receiver = as.character(em$receiver), increment = 1)
mail <- linkEvents(x = mail, changeEvent = events, nodes = nodes)
dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = mail)
# (a) event network, all ties created at time 0
hev <- data.frame(time = 0, sender = as.character(hi$sender), receiver = as.character(hi$receiver), increment = 1)
hier <- defineNetwork(nodes = nodes, directed = TRUE)
hier <- linkEvents(x = hier, changeEvent = hev, nodes = nodes)
ma <- estimate(dep ~ inertia(mail) + recip(mail) + indeg(mail) + tie(hier), model = "DyNAM", subModel = "choice")
# (b) static network
M <- matrix(0, length(ea), length(ea), dimnames = list(ea, ea))
for (k in seq_len(nrow(hi))) M[as.character(hi$sender[k]), as.character(hi$receiver[k])] <- 1
hs <- defineNetwork(M, nodes = nodes, directed = TRUE)
mb <- estimate(dep ~ inertia(mail) + recip(mail) + indeg(mail) + tie(hs), model = "DyNAM", subModel = "choice")
for (m in list(ma, mb)) print(c(coef(m), ll = as.numeric(summary(m)$logLikelihood)), digits = 15)
