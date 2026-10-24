# dev/dynam_unit21c_rateintercept_weighted_crosscheck.R -- WITH-INTERCEPT
# rate sub-model combined with weighted=TRUE (batch 18). Real timestamps,
# same toy network shape as unit21b.
library(goldfish)
set.seed(13)
n <- 6
actors <- data.frame(label = paste0("A", 1:n), present = TRUE)
nodes <- defineNodes(actors)
mat0 <- matrix(0, n, n, dimnames = list(actors$label, actors$label))
diag(mat0) <- NA
net <- defineNetwork(matrix = mat0, nodes = nodes, directed = TRUE)

nevents <- 14
events <- data.frame(
  time = cumsum(round(runif(nevents, 1, 10), 2)),
  sender = paste0("A", sample(1:n, nevents, replace = TRUE)),
  receiver = paste0("A", sample(1:n, nevents, replace = TRUE))
)
events <- events[events$sender != events$receiver, ]
events$increment <- 1
net <- linkEvents(x = net, changeEvents = events, nodes = nodes)
depnet <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)

mod <- estimate(depnet ~ 1 + outdeg(net, weighted = TRUE), model = "DyNAM", subModel = "rate")
print(summary(mod))
cat(sprintf("INTERCEPT=%.8f OUTDEGW=%.8f LOGLIK=%.8f\n",
    coef(mod)["Intercept"], coef(mod)["outdeg"], as.numeric(logLik(mod))))

write.csv(events, "dev/dynam_crosscheck_rateintweighted_events.csv", row.names = FALSE)
