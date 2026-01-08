MM1_sim <- function(lambda, mu, T_end){
  Z <- c()
  n <- c()
  Z_type <- c()
  Z[1] <- rexp(1, rate=lambda)
  Z_type[1] <- 1
  n[1] <- 272
  i=1
  while (sum(Z) < T_end){
    i=i+1
    # Either birth or death can occur if the population size is not 0
    if (rev(n)[1] > 0){
      Z[i] <- rexp(1, rate=lambda+mu)
      Z_type[i] <- sample(c(-1,1), size=1, prob=c(mu/(lambda+mu), lambda/(lambda+mu))) 
      n[i] <- n[i-1] + Z_type[i]
    }
    # only birth can occur if the population size is 0
    else if (rev(n)[1] == 0){
      Z[i] <- rexp(1, rate=lambda)
      Z_type[i] <- 1
      n[i] <- n[i-1] + Z_type[i]
    }
  }
  return(list(Z=Z, n=n, Z_type=Z_type))
  
}

set.seed(123)
par(mfrow=c(2,2))
lambda0 = c(0)
mu0 = c(0.26)
for (i in 1:4){
  res <- MM1_sim(lambda=lambda0[i], mu=mu0[i], T_end=19)
  plot(cumsum(res$Z), res$n, type='l', xlab='time', ylab='population', main=paste0('lambda=', lambda0[i], ',', 'mu=', mu0[i]))
}