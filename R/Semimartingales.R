#1.2) Continuous Semimartingales
#1.2.1) Continuous Martingales

# Brownian motion
BB <- function(seed, delta, tau){
  set.seed(seed)
  time <- seq(from = 0, to = tau, by = delta)
  nticks <- length(time) - 1 # number of sample points
  std <- sqrt(delta) # standard deviation Brownian motion
  b <- cumsum (c(0, rnorm (nticks , sd = std))) # b is a simulated Brownian motion path sampled at the tick marks.
  b
}

# Compensated Squared Brownian motion
CBB <- function(seed ,delta ,tau){
  time <- seq(from = 0, to = tau, by = delta)
  comb <- BB(seed, delta, tau) ^ 2 - time
  comb
}

# Exponential Brownian motion
EBB <- function(seed, delta, tau){
  time <- seq(from = 0, to = tau, by = delta)
  ex <- exp(BB(seed, delta, tau) - 1 / 2 * time)
  ex
}

#1.2.2) Continuous non - martingale semimartingales
# Brownian motion with drift
DBB <- function(seed, drift, delta, tau){
  time <- seq (from = 0, to = tau, by = delta)
  dbb <- BB(seed, delta, tau ) + drift * time
  dbb
}

# Geom.Brownian motion with start in 1 (0 is degenerated case )
GBB <- function(seed, mu, sigma, delta, tau){
  time <- seq(from = 0, to = tau, by = delta)
  gbb <- 1 * exp((mu - sigma ^ 2 / 2) * time + sigma * BB(seed, delta, tau))
  gbb
}

# Ornstein Uhlenbeck
ornstein_uhlenbeck <- function(T, n, nu, lambda, sigma, x0){
  dw <- rnorm(n, 0, sqrt(T / n))
  dt <- T / n
  x <- c(x0)
  for (i in 2:(n + 1)){
    x[i] <- x[i-1] + lambda * (nu - x[i-1]) * dt + sigma * dw[-1]
  }
  return (x);
}  
  