#1.) Continuous testing
#1.1) Functions

# Sample quadratic variation ( realised volatility )
fqv <- function(process){
  process <- as.numeric(process)
  y <- cumsum(diff(process)^2)
  c(y, y[length(y)])
}

# New grid of new mesh DELTA on [0, fqv ( process )[ length ( process )]]
gitterNeu <- function(process, Delta){
  seq(from = 0, to = fqv(process)[length(process)], by = Delta)
}

# Time changed process ( approximated version like in Wittmann )
yNeuV <- function(process, gitterNeu){
  N <- length(gitterNeu)
  idx <- sapply(1:N, function(i){
    which.min(abs(fqv(process) - i * gitterNeu[2])) # gitterNeu [2] == Delta
   })
  process[idx]
}

# Optimal Delta based on smallest ks. statistic
optimalDelta <- function(process, rate, upper){
  interval <- seq(from = 0 + rate, to = upper, by = rate) #We choose rate =\ delta /4 and upper =1/2
  x <- c()
  for (i in 1:length(interval)) {
    x[i] <- as.vector(ks.test(1 / sqrt(interval[i]) * diff(yNeuV(process, gitterNeu(process, interval[i]))), pnorm)$statistic)
    if (length(1 / sqrt(interval[i]) * diff(yNeuV(process, gitterNeu(process, interval[i])))) == 1) break # Avoid error if N=1
  }
 interval[which.min(x)]
}
