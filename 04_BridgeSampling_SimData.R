#-------------------------------Calculating Bayes' Factor with Bridge Sampling--------------------------------------------------

# Author: Hamburger Deern
# Date: 24.06.2020

# Goal: Calculating the Bayes' Factors with Bridge Sampling

# Package Documentation: https://cran.r-project.org/web/packages/bridgesampling/bridgesampling.pdf

#-----------------------------Preliminaries---------------------------------------------------

library(beepr)
library(R2jags)
library(bridgesampling)
library(LaplacesDemon)

#-----------------------------Loading the data-----------------------------------------------

#-------------------------Option 2: Simulated data--------------------------------------------

# Put in the path of the csv-file with the simulated data:
simData <- read.csv("00_Data_Simulation.csv", row.names = 1)

# Select the distribution by changing the accessed column:
y <- simData$sim_genpois

#-----------------------------Fitting the models---------------------------------------------

# Set the hypothesis:
H0 <- "GenPois"
H1 <- "NegBinom"

set.seed(1)
getSamplesModel_Pois <- function(data) {
  
  model <- "
    model{
      for (i in 1:N) {
        y[i] ~ dpois(lambda)
      }
      lambda ~ dgamma(1,1)
    }"
  
  s <- jags(data, parameters.to.save = c("lambda"),
            model.file = textConnection(model),
            n.chains = 1, n.iter = 1000, n.burnin = 100)
  
  return(s)
}

getSamplesModel_NegBinom <- function(data) {
  
  model <- "
    model{
      for (i in 1:N) {
        y[i] ~ dnegbin(p[i],theta)
        p[i] <- theta/(theta+lambda)
      }
      theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
      omega ~ dbeta(1,1)
      lambda ~ dgamma(1,1)
    }"
  
  s <- jags(data, parameters.to.save = c("lambda", "theta"),
            model.file = textConnection(model),
            n.chains = 1, n.iter = 1000, n.burnin = 100)
  
  return(s)
}

getSamplesModel_GenPois <- function(data) {
  
  model <- "
    data {
      C <- 10000
      for (i in 1:N) {
        ones[i] <- 1
      } 
    }
    model {
     for (i in 1:N) {
        spy[i] <- (((1-omega)*lambda*((1-omega)*lambda + omega*y[i])^(y[i]-1)*exp(-1*((1-omega)*lambda + omega*y[i]))) / exp(logfact(y[i]))) / C
        ones[i] ~ dbern(spy[i])
      }
    omega ~ dbeta(1,1)
    lambda ~ dgamma(1,1)
  }"   
  
  s <- jags(data, parameters.to.save = c("lambda","omega"),
            model.file = textConnection(model),
            n.chains = 1, n.iter = 1000, n.burnin = 100)
  
  return(s)
}

#---------------------------Specifying the Unnormalized Log Posterior Function----------------------------
  
### get posterior samples ###
# create data list for JAGS
data <- list(y = y, N = length(y))

# fit models
if (H0 == "Pois"){
    samples_H0 <- getSamplesModel_Pois(data)
} else if (H0 == "NegBinom"){   
    samples_H0 <- getSamplesModel_NegBinom(data)
} else if (H0 == "GenPois"){  
    samples_H0 <- getSamplesModel_GenPois(data)
}

if (H1 == "Pois"){
  samples_H1 <- getSamplesModel_Pois(data)
} else if (H1 == "NegBinom"){   
  samples_H1 <- getSamplesModel_NegBinom(data)
} else if (H1 == "GenPois"){  
  samples_H1 <- getSamplesModel_GenPois(data)
}

### functions for evaluating the unnormalized posteriors on log scale ###

log_posterior_Pois <- function(samples.row, data) {
  
  lambda <- samples.row[[ "lambda" ]]
  
  sum(dpois(data$y, lambda, log = TRUE))
}

log_posterior_NegBinom <- function(samples.row, data) {
  
  lambda <- samples.row[[ "lambda" ]]
  theta <- samples.row[[ "theta" ]]
  
  sum(dnbinom(data$y, theta, prob = theta/(theta+lambda), log = TRUE))
}

log_posterior_GenPois <- function(samples.row, data) {
  
  lambda <- samples.row[[ "lambda" ]]
  omega <- samples.row[[ "omega" ]]
  
  sum(dgpois(data$y, lambda, omega, log = TRUE))
}

if (H0 == "Pois"){
  log_posterior_H0 <- log_posterior_Pois
} else if (H0 == "NegBinom"){   
  log_posterior_H0 <- log_posterior_NegBinom
} else if (H0 == "GenPois"){  
  log_posterior_H0 <- log_posterior_GenPois
}

if (H1 == "Pois"){
  log_posterior_H1 <- log_posterior_Pois
} else if (H1 == "NegBinom"){   
  log_posterior_H1 <- log_posterior_NegBinom
} else if (H1 == "GenPois"){  
  log_posterior_H1 <- log_posterior_GenPois
}

#---------------------------------Specifying the Parameter Bounds--------------------------------------

# specify parameter bounds H0
cn <- colnames(samples_H0$BUGSoutput$sims.matrix)
cn <- cn[cn != "deviance"]
lb_H0 <- rep(-Inf, length(cn))
ub_H0 <- rep(Inf, length(cn))
names(lb_H0) <- names(ub_H0) <- cn
lb_H0[[ "lambda" ]] <- 0

if (H0 == "NegBinom"){   
  lb_H0[[ "theta" ]] <- 0
} else if (H0 == "GenPois"){  
  lb_H0[[ "theta" ]] <- 0
}

# specify parameter bounds
cn <- colnames(samples_H1$BUGSoutput$sims.matrix)
cn <- cn[cn != "deviance"]
lb_H1 <- rep(-Inf, length(cn))
ub_H1 <- rep(Inf, length(cn))
names(lb_H1) <- names(ub_H1) <- cn
lb_H1[[ "lambda" ]] <- 0

if (H1 == "NegBinom"){   
  lb_H1[[ "theta" ]] <- 0
} else if (H1 == "GenPois"){  
  lb_H1[[ "theta" ]] <- 0
}

#---------------------------------Computing the (Log) Marginal Likelihoods--------------------------

# compute log marginal likelihood via bridge sampling for H0
H0.bridge <- bridge_sampler(samples = samples_H0, data = data,
                            log_posterior = log_posterior_H0, lb = lb_H0,
                            ub = ub_H0, silent = TRUE)
# compute log marginal likelihood via bridge sampling for H1
H1.bridge <- bridge_sampler(samples = samples_H1, data = data,
                            log_posterior = log_posterior_H1, lb = lb_H1,
                            ub = ub_H1, silent = TRUE)

#------------------------------------Bayesian Model Comparison-----------------------------------

# compute Bayes factor
# LOG
BF01 <- bf(H0.bridge, H1.bridge)
print(log10(BF01$bf))
BF02 <- bf(H1.bridge, H0.bridge)
print(log10(BF02$bf))

# LN
BF01 <- bf(H0.bridge, H1.bridge, log = TRUE)
print(BF01)
BF02 <- bf(H1.bridge, H0.bridge, log = TRUE)
print(BF02)

# DIC
samples_H0$BUGSoutput$DIC
samples_H1$BUGSoutput$DIC

beep(1)

# Error
error_measures(H0.bridge, H1.bridge)
error_measures(H1.bridge, H1.bridge)
