#2.3) Asymptotic Distribution
# Simulating the asymptotic test statistic of S_n and S_n * like in Park and
# Whang (2005) and Phillips and Jin (2014)
# -->Ito - sum approach

asymdistr.PW <- function(NSIMS, boundaries, rate){
  set.seed(1)
  delta <- 1 / 1000
  time <- seq(from = 0, to = 1, by = delta)
  nticks <- length(time) - 1 # number of sample points
  std <- sqrt(delta) # standard deviation Brownian motion
  x <- seq(from = -boundaries, to = boundaries, by = rate) # Simulations show that boundaries =5 is appropriate . We choose rate =0.01
  M <- c()
  S <- c()
  Ergebnisse <- list()
  for (i in 1:NSIMS) {
    b <- cumsum(c(0, rnorm(nticks, sd = std))) # b is a simulated Brownian motion path sampled at the tick marks .
    for (j in 1:length(x)) {
      ito.integral <- 0 # simulating the ito - integral 1*( B(s) <=x) dBs from 0 to  1
      for (k in 1:(length(b) - 1)) {
        ito.integral <- ito.integral + 1 * (b[k] <= x[j]) * (b[k + 1] - b[k])
      }
      M[j] <- ito.integral
    }
    Ergebnisse[[1]] <- M
    S[i] <- max(abs(M))
  }
  Ergebnisse[[2]] <- S
  Ergebnisse
}

# Simulating the asymptotic teststatistic of the martingale - test of ( modified )
#Park and Whang with a Karhunen -Loeve - approach
asymdistrKL.PW <- function(){
  NSIMS <- 1000 # number of simulations
  T <- 1 # number of sample points between 0 and 1
  N <- 1000
  t <- seq(0, T, length = N + 1)
  b <- numeric(N + 1)
  x <- seq(from = -5, to = 5, by = 0.01)
  Z <- rnorm(NSIMS)
  phi <- function(i, t, T){
    (2 * sqrt(2 * T)) / ((2 * i + 1) * pi) * sin(((2 * i + 1) * pi * t) / (2 * T)) # basis of orthogonal functions in L2 ([0 ,T])
  }
  S <- c()
  M <- c()
  for (i in 1:NSIMS) {
    for (j in 1:length(x)) {
      for (i in 2:(N + 1)) {
        b[i] <- sum(Z * sapply(1:NSIMS, function(x) phi(x, t[i], T)))
      }# b is a simulated Brownian motion path , based on Karhunen - Love approach
      ito.integral <- 0 # simulating the ito - integral 1*( B(s) <=x) dBs from 0 to 1
      for (k in 1:(length(b) - 1)) {
        ito.integral <- ito.integral + 1 * (b[k] <= x[j]) * (b[k + 1] - b[k])
      }
      M[j] <- ito.integral
    }
    M
    S[i] <- max(abs(M)) # supremum -> generalized KS - test
  }
  S
}