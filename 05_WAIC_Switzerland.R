# ---- Probabilistic Modelling ---------------------------------------------

# Author: Hamburger Deern

# Goal: Calculate WAIC for the Switzerland data
## M1: Poisson
## M2: Negative binomial
## M3: Generalized poisson

# Documentation https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/

# ---- 0. Preliminaries ---------------------------------------------------

library(lavaan)
library(loo)
library(R2jags)
library(LaplacesDemon)

set.seed(1)

# Import Data
claim_data <- read.csv("claim_data.csv", row.names = 1)

y <- claim_data$Switzerland_1961

y_long <- c()
for (i in 1:length(y)){
  claims = i-1
  y_long <- c(y_long, rep(claims, times = y[i]))
}
y_long <- as.integer(y_long)

dataList <- list("y" = y_long, "N" = length(y_long))

# Run JAGS model ----
n.burnin = 1000
n.iter = 2000
n.chains = 2

# Model 1: Poisson
cat("model{
  for (i in 1:N) {
    y[i] ~ dpois(lambda)
    #loglik[i] <- logdensity.pois(y[i], lambda) # add log likelihood computation for each observation!
    loglik[i] <- log(dpois(y[i], lambda))
  }
  lambda ~ dgamma(0.0001,0.0001)
}", file='claim_model.txt')

s_m1 <- jags(data = dataList, #5.4 GB
          parameters.to.save = c("lambda", "loglik"),
          model.file = 'claim_model.txt',
          n.chains = n.chains, 
          n.iter = n.iter, 
          n.burnin = n.burnin)

loglik_m1 <- s_m1$BUGSoutput$sims.list$loglik #1.8 GB
size_m1 <- object.size(s_m1)
rm(s_m1) # remove JAGS output to save memory

# Model 2: Negative Binomial
cat("model{
      for (i in 1:N) {
        y[i] ~ dnegbin(p[i],theta)
        p[i] <- theta/(theta+lambda)
        loglik[i] <- log(dnegbin(y[i], p[i], theta))

      }
      theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
      omega ~ dbeta(1,1)
      lambda ~ dgamma(0.0001,0.0001)
}", file='claim_model.txt')

s_m2 <- jags(data = dataList, #5.4
             parameters.to.save = c("lambda", "theta", "loglik"),
             model.file = 'claim_model.txt',
             n.chains = n.chains, 
             n.iter = n.iter, 
             n.burnin = n.burnin)

loglik_m2 <- s_m2$BUGSoutput$sims.list$loglik
size_m2 <- object.size(s_m2)
rm(s_m2) # remove JAGS output to save memory

# Model 3: Generalized Poisson
cat("data {
  C <- 10000
  for (i in 1:N) {
    ones[i] <- 1
  }
} model {
   for (i in 1:N) {
      spy[i] <- (((1-omega)*lambda*((1-omega)*lambda + omega*y[i])^(y[i]-1)*exp(-1*((1-omega)*lambda + omega*y[i]))) / exp(logfact(y[i]))) / C
      ones[i] ~ dbern(spy[i])
      loglik[i] <- log((((1-omega)*lambda*((1-omega)*lambda + omega*y[i])^(y[i]-1)*exp(-1*((1-omega)*lambda + omega*y[i]))) / exp(logfact(y[i]))))
   }
  omega ~ dbeta(1,1)
  lambda ~ dgamma(0.0001,0.0001)
}", file='claim_model.txt')


s_m3 <- jags(data = dataList, #5.4
             parameters.to.save = c("lambda", "omega", "loglik"),
             model.file = 'claim_model.txt',
             n.chains = n.chains, 
             n.iter = n.iter, 
             n.burnin = n.burnin)

loglik_m3 <- s_m3$BUGSoutput$sims.list$loglik
size_m3 <- object.size(s_m3)
rm(s_m3) # remove JAGS output to save memory

# Calculate WAIC
waic_m1 <- waic(loglik_m1)
waic_m2 <- waic(loglik_m2)
waic_m3 <- waic(loglik_m3)

loo_compare(waic_m1, waic_m2, waic_m3)

Sys.time()
