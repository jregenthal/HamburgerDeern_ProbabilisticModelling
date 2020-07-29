# ---- Probabilistic Modelling ---------------------------------------------

# Author: Hamburger Deern

# Goal: Post-Processing RJMCMC for Poisson Data
## M1: Poisson
## M2: Negative binomial
## M3: Generalized poisson

# Package Documentation: https://cran.r-project.org/web/packages/rjmcmc/rjmcmc.pdf

# ---- 0. Preliminaries ---------------------------------------------------
library(madness)
library(rjmcmc)
library(LaplacesDemon)
library(ggplot2)

# Get the function diagMCMC for diagnositics
source("DBDA2E-utilities.R")
set.seed(1)

# ---- 1. Import data ----------------------------------------------

claim_data <- read.csv("00_Data_Simulation.csv", row.names = 1)

y_long <- claim_data$sim_pois


# ---- 2. Iterations  ---------------------------------------
n.adapt = 500
n.burn = 1000
n.iter = 10000
n.chains = 3
n.rjmcmc = 10000


## ---- 2a. M1: Poisson model ---------------------------------------------
# Define the model
modelString_pois <- "model{
  for (i in 1:N) {
    y[i] ~ dpois(lambda)
  }
  lambda ~ dgamma(1,1)
}"

writeLines(modelString_pois, con = "claim_model.txt")

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_pois <- jags.model(file = "claim_model.txt",
                             data = dataList,
                             n.chains = n.chains, 
                             n.adapt = n.adapt)

update(jagsModel_pois, n.burn)

# Simulate the model
codaSamples_pois = coda.samples(jagsModel_pois, 
                                variable.names = c("lambda"), 
                                n.iter = n.iter)
summary(codaSamples_pois)
m1_lambda = mean(codaSamples_pois[[1]][,1])

chain_poisson <- data.frame(iter = 1:n.iter, codaSamples_pois[[1]])


## ---- 2b. M2: Negative binomial model ---------------------------------------------
# Define the model
modelString_negbin <- "model{
    for (i in 1:N) {
      y[i] ~ dnegbin(p[i],theta)
      p[i] <- theta/(theta+lambda)
    }
    theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
    omega ~ dbeta(1,1)
    lambda ~ dgamma(1,1)
}"

writeLines(modelString_negbin, con = "claim_model.txt")

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_negbin <- jags.model(file = "claim_model.txt",
                               data = dataList,
                               n.chains = n.chains, 
                               n.adapt = n.adapt)

update(jagsModel_negbin, n.burn)

# Simulate the model
codaSamples_negbin = coda.samples(jagsModel_negbin, 
                                  variable.names = c("lambda", "theta"), 
                                  n.iter = n.iter)
summary(codaSamples_negbin)
m2_lambda = mean(codaSamples_negbin[[1]][,1])
m2_theta = mean(codaSamples_negbin[[1]][,2])

chain_negbin <- data.frame(iter = 1:n.iter, codaSamples_negbin[[1]])


## ---- 2c. M3: Generalized poisson model ---------------------------------------------
# Define the model
modelString_genpois <- "data {
  C <- 10000
  for (i in 1:N) {
    ones[i] <- 1
  }
} model {
   for (i in 1:N) {
      spy[i] <- (((1-omega)*lambda*((1-omega)*lambda + omega*y[i])^(y[i]-1)*exp(-1*((1-omega)*lambda + omega*y[i]))) / exp(logfact(y[i]))) / C
      ones[i] ~ dbern(spy[i])
    }
  omega ~ dbeta(1,1)
  lambda ~ dgamma(1,1)
}"

writeLines(modelString_genpois, con = "claim_model.txt")

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_genpois <- jags.model(file = "claim_model.txt",
                               data = dataList,
                               n.chains = n.chains, 
                               n.adapt = n.adapt)

update(jagsModel_genpois, n.burn)

# Simulate the model
codaSamples_genpois = coda.samples(jagsModel_genpois, 
                                  variable.names = c("lambda", "omega"), 
                                  n.iter = n.iter)
summary(codaSamples_genpois)
m3_lambda = mean(codaSamples_genpois[[1]][,1])
m3_omega = mean(codaSamples_genpois[[1]][,2])

chain_genPois <- data.frame(iter = 1:n.iter, codaSamples_genpois[[1]])


# ---- 3 RJMCMC pre-definitions --------------------------------------
# M1: Poisson
## Likelihood
L1=function(parameter){sum(dpois(y_long,
                                 parameter[1], 
                                 log = FALSE))}
## Prior
p.prior1=function(parameter){
  dgamma(parameter[1], 0.0001, 0.0001,log=FALSE)} 

# M2: Negative binomial
## https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/NegBinomial
## x = y, size = theta, prob = theta / (theta+lambda)
## Likelihood
L2=function(parameter){
  return(sum(dnbinom(y_long,
                     parameter[2],
                     prob = parameter[2]/(parameter[2]+parameter[1]),
                     log = FALSE)))}
## Prior
p.prior2=function(parameter){
  dgamma(parameter[1], 1, 1,log=FALSE)
  + 1/2 * parameter[1] * parameter[2]^(-2) * (1 + parameter[1]/parameter[2])^(-1.5)} 

# M3: Generalized Poisson
## Likelihood
L3=function(parameter){
  return(sum(dgpois(y_long,
                    parameter[1],
                    parameter[2],
                    log = FALSE)))}
## Prior
p.prior3=function(parameter){
  dgamma(parameter[1], 1, 1,log=FALSE)
  + dbeta(parameter[2], 1 ,1, log=FALSE)
} 


# ---- 5. Draw from posteriors -----------------------------
# Sample from posterior

C1 <- as.matrix(chain_poisson[,2])
C2 <- as.matrix(chain_negbin[,2:3])
C3 <- as.matrix(chain_genPois[,2:3])

# ---- 6. Run post-processing RJMCMC --------------------------------------
# M2 is left out as it only throws errors as theta does not have a realistic value
out = defaultpost(posterior = list(C1, C3), # JAGS MCMC output as matrices
                  likelihood = list(L1, L3), # Functions to sample from model likelihoods
                  param.prior = list(p.prior1, p.prior3), # Functions to sample from priors
                  model.prior = c(0.5, 0.5), # Model probabilities
                  chainlength = n.rjmcmc,
                  save.all = TRUE) 

summary(out)

# ln Bayes factors
logNAT_BayesFactor = log(out[[1]]$`Bayes Factors`)
print(logNAT_BayesFactor)

# log10 Bayes factors
log10_BayesFactor = log10(out[[1]]$`Bayes Factors`)
print(log10_BayesFactor)


# ---- 7. Visualization --------------------------------
# Plot of posterior model probabilities
plot(out)

# Examine Priors, Likelihoods, Posteriors
densities <- data.frame(iter=1:n.rjmcmc, out$densities[[1]])

ggplot(data = densities) + 
  geom_line(aes(x = iter, y = Posterior.M1, color = 'red'), show.legend = T) +
  geom_line(aes(x = iter, y = Posterior.M2, color = 'blue'), show.legend = T) +
  scale_color_discrete(name = "Models", labels = c("M1: Poisson", "M3: Generalized Poisson"))

ggplot(data = densities) + 
  geom_line(aes(x = iter, y = Likelihood.M1, color = 'red'), show.legend = T) +
  geom_line(aes(x = iter, y = Likelihood.M2, color = 'blue'), show.legend = T) +
  scale_color_discrete(name = "Models", labels = c("M1: Poisson", "M3: Generalized Poisson"))

ggplot(data = densities) + 
  geom_line(aes(x = iter, y = Prior.M1, color = 'red'), show.legend = T) +
  geom_line(aes(x = iter, y = Prior.M2, color = 'blue'), show.legend = T) +
  scale_color_discrete(name = "Models", labels = c("M1: Poisson", "M3: Generalized Poisson"))