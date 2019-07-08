#1.3) Evaluating size performance and power performance . Example of Brownian motion
#1.3.1) Size performance , based on the KS - test like in Vasudev (2007)

powerBB <- function(rate, nsim, alpha, tau, delta ,upper){
  p <- c()
  for (i in 1:nsim) {
    time <- seq(from = 0, to = tau, by = delta)
    nticks <- length(time) - 1 # number of sample points
    std <- sqrt(delta) # standard deviation of Brownian motion
    b <- cumsum(c(0, rnorm(nticks, sd = std))) #b is Brownian motion path
    x <- optimalDelta(b, rate, upper)
    p[i] <- as.vector(ks.test(1 / sqrt(x) * diff(yNeuV(b, gitterNeu(b, x))), pnorm)$p.value)
  }
  p
  pow <- c()
  for (i in 1:length(p)) {
    if (p[i] < alpha) {
      pow[i] <- 1 
    } else {
        pow[i] <- 0
        }
    }
    sum(pow) / length(p)
}

#1.3.2) Size performance , based on the HZ - approach
powerBBHZ <- function(rate, nsim, alpha, tau, delta, upper){
  p <- c()
  for (i in 1:nsim) {
    time <- seq(from = 0, to = tau, by = delta)
    nticks <- length(time) - 1 # number of sample points
    std <- sqrt(delta) # standard deviation of Brownian motion
    b <- cumsum(c(0, rnorm(nticks, sd = std))) #b is Brownian motion path
    od <- optimalDelta(b, rate, upper) # Note : Optimal Delta is chosen based on KS - test statistic (!)
    y <- 1 / sqrt(od) * diff(yNeuV(b, gitterNeu(b, od)))
    if (length(y) %% 2 != 0) {
      y <- y[-length(y)] # obtain two vectors of the same length
    }
    odd_indexes <- seq(1, length(y), 2) # split sample to obtain two vectors of standardized increments
    even_indexes <- seq(2, length(y), 2)
    p[i] <- as.vector(hzTest(cbind(y[odd_indexes], y[even_indexes]))@p.value)
  }
  p
  pow <- c()
  for (i in 1:length(p)) {
    if (p[i] < alpha) {
      pow[i] <- 1
    } else {
        pow[i] <- 0
        }
    }
  sum(pow) / length(p)
}