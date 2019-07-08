#2.) Discrete Testing
#2.1) Functions

#Qn(x) like in Park and Whang (2005) , plot version
qn <- function(process, mesh){
  Qn <- c()
  sigma <- sqrt(1 / (length(process) - 1) * sum(diff(process) ^ 2)) # process_0 =0
  x <- seq(from = min(process / sigma) - mesh, to = max(process / sigma ) + mesh, by = mesh)
    for (i in 1:length(x)) {
    Qn[i] <- 1 / (sqrt(length(process) - 1)) * (sum((diff(process) / sigma) * 1 * (process[-length(process)] / sigma <= x[i])))
    }
 plot(Qn, cex = .4, type = "l")
}

# Test statistic S_n like in Park and Whang (2005)
Sn <- function(process, mesh){
  Qn <- c()
  Sn <- c()
  sigma <- sqrt(1 / (length(process) - 1) * sum(diff(process) ^ 2) ) # process_0 =0
  x <- seq( from = min( process / sigma ) - mesh , to = max(process/sigma) + mesh, by = mesh)
  for (i in 1:length(x)) {
    Qn[i] <- 1 / (sqrt(length(process) - 1)) * (sum((diff(process) / sigma) * 1 * (process[-length(process)] / sigma <= x[i])))
  }
  max(abs(Qn))
}

# Test statistic S_n * based on modified Park and Whang (2005) test like in
#Phillips and Jin (2014)
S_n_star <- function(process, mesh){
  reg <- lm(process[-1] ~ process[-length(process)]) # regress Y_t on Y_t -1 <-> eliminate Y_0 =0 resp . the last value
  ut <- reg$residuals # standardization based on OLS residual
  # sigmahead <- sqrt (1 / (length(process) - 1) * sum(ut ^ 2))
  Jn <- c()
  x <- seq(from = min(process) - mesh, to = max(process) + mesh, by = mesh)
  for (i in 1:length(x)) {
    Jn[i] <- (sum(diff(process) * 1 * (process[-length(process)] <= x[i]))) / (sqrt(sum(ut ^ 2)))
  }
  max(abs(Jn))
}

# Test statistic based on modified Park and Whang (2005) test with drift like in
#Phillips and Jin (2014) , not discussed in thesis
GKSstar <- function(process, mesh){
  reg <- lm(process[-1] ~ process[-length(process)])
  ut <- reg$residuals
  # sigmahead <- sqrt(1 / (length(process) - 1) * sum(ut ^ 2))
  JnS <- c()
  x <- seq(from = min(process) - mesh, to = max(process) + mesh, by = mesh)
  for (i in 1:length(x)) {
    JnS[i] <- (sum((diff(process) - 1 / (length(process) - 1) * sum(diff(process))) * 1 * (process[-length(process)] <= x[i]))) / (sqrt(sum(ut ^ 2)))
  }
  max(abs(JnS))
}

#2.2) Simulating processes under the null and under the alternative , possible with GARCH - errors
random.walk.GARCH <- function(seed, n, arch, garch){
  n <- n + 1
  set.seed(seed)
  e <- rnorm(n, 0, 1)
  u <- numeric(); sigsq <- numeric(); x <- numeric()
  u[1] <- 0
  sigsq[1] <- 0
  x[1] <- 0
  for (i in 2:n) {
    sigsq[i] <- 1 + arch * (u[i - 1] ^ 2) + garch * (sigsq[i - 1])
    u[i] <- sqrt(sigsq[i]) *  e[i]
    x[i] <- x[i - 1] + u[i]
  }
  x
}

#Simulating ARMA (1 ,1) process
arma <- function(seed, n, ar, ma){
  n <- n + 1
  set.seed(seed)
  e <- rnorm(n, 0, 1)
  x <- numeric()
  x[1] <- 0
  e[1] <- 0
  for (i in 2:n) {
    x[i] <- ar * x[i - 1] + ma * e[i - 1] + e[i]
  }
  x
}