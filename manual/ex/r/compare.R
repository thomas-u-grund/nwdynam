# goldfish side of the comparison of chapter 11: fits the models on the two event logs with goldfish and writes
# ex/r/goldfish_results.csv (model, term, estimate, se, loglik). Run from manual/ex/r. Needs R with goldfish 1.6.12.
suppressMessages(library(goldfish))
fit <- function(tag, ev, actors, effects, sub) {
  nodes <- defineNodes(data.frame(label = actors))
  net <- defineNetwork(nodes = nodes, directed = TRUE)
  events <- data.frame(time = ev$t, sender = as.character(ev$sender), receiver = as.character(ev$receiver), increment = 1)
  net <- linkEvents(x = net, changeEvent = events, nodes = nodes)
  dep <- defineDependentEvents(events = events, nodes = nodes, defaultNetwork = net)
  assign("dep", dep, envir = globalenv()); assign("net", net, envir = globalenv()); assign("nodes", nodes, envir = globalenv()); assign("events", events, envir = globalenv())
  f <- as.formula(paste("dep ~", effects), env = globalenv())
  m <- estimate(f, model = "DyNAM", subModel = sub)
  s <- summary(m)
  data.frame(model = tag, term = names(coef(m)), estimate = as.numeric(coef(m)), se = as.numeric(s$standardErrors), loglik = as.numeric(s$logLikelihood))
}
em <- read.csv("emails.csv"); ea <- read.csv("email_actors.csv")$label
so <- read.csv("socevo.csv"); sa <- as.character(1:84)
res <- rbind(
  fit("soc_choice", so, sa, "inertia(net) + recip(net) + indeg(net)", "choice"),
  fit("soc_rate", so, sa, "indeg(net) + outdeg(net)", "rate"),
  fit("soc_rate_int", so, sa, "1 + indeg(net) + outdeg(net)", "rate"),
  fit("em_choice", em, ea, "inertia(net) + recip(net) + indeg(net)", "choice"),
  fit("em_choice_trans", em, ea, "inertia(net) + recip(net) + indeg(net) + trans(net)", "choice"),
  fit("em_choice_win", em, ea, "inertia(net, window = 0.5) + recip(net) + indeg(net)", "choice"),
  fit("em_rate", em, ea, "indeg(net) + outdeg(net)", "rate"),
  fit("em_rate_int", em, ea, "1 + indeg(net) + outdeg(net)", "rate")
)
write.csv(res, "goldfish_results.csv", row.names = FALSE)
print(res)
