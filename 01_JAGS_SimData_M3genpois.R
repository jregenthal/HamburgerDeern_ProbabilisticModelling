# ---- Probabilistic Modelling ---------------------------------------------

# Author: Hamburger Deern
# Date: 28.06.2020

# Goal: Run MCMC for three models
## M1: Poisson
## M2: Negative binomial
## M3: Generalized poisson

# ---- 0. Preliminaries ---------------------------------------------------
library(LaplacesDemon)
library(ggplot2)

# Get the function diagMCMC for diagnositics
source("DBDA2E-utilities.R")
set.seed(1)

# ---- 1. Import data ----------------------------------------------

claim_data <- read.csv("00_Data_Simulation.csv", row.names = 1)

y_long <- claim_data$sim_genpois

# ---- 2. Run MCMC models ---------------------------------------
n.adapt = 500
n.burn = 1000
n.iter = 10000
n.chains = 3


## ---- 2a. M1: Poisson model ---------------------------------------------

# Define the model
modelString_pois <- "model{
  for (i in 1:N) {
    y[i] ~ dpois(lambda)
  }
  lambda ~ dgamma(1,1)
}"

writeLines(modelString_pois, con = path_MCMCtxt)

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_pois <- jags.model(file = path_MCMCtxt,
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

# Visualizations

# Convergence diagnostics:
diagMCMC( codaObject=codaSamples_pois, parName="lambda" )

# Posterior descriptives:
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples_pois[,"lambda"] , main="lambda" , xlab=bquote(lambda) )

# Re-plot with different annotations:
plotPost( codaSamples_pois[,"lambda"] , main="lambda" , xlab=bquote(lambda) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )

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

writeLines(modelString_negbin, con = path_MCMCtxt)

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_negbin <- jags.model(file = path_MCMCtxt,
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

# Visualization

# Convergence diagnostics:
diagMCMC( codaObject=codaSamples_negbin, parName="lambda" )

# Posterior descriptives:
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples_negbin[,"lambda"] , main="lambda" , xlab=bquote(lambda) )

# Re-plot with different annotations:
plotPost( codaSamples_negbin[,"lambda"] , main="lambda" , xlab=bquote(lambda) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )

# Convergence diagnostics:
diagMCMC(codaObject=codaSamples_negbin, parName="theta" )

# Posterior descriptives:
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples_negbin[,"theta"] , main="theta" , xlab=bquote(theta) )

# Re-plot with different annotations:
plotPost( codaSamples_negbin[,"theta"] , main="theta" , xlab=bquote(theta) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )

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

writeLines(modelString_genpois, con = path_MCMCtxt)

# Compile the model
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel_genpois <- jags.model(file = path_MCMCtxt,
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

# Visualizations

# Convergence diagnostics:
diagMCMC( codaObject=codaSamples_genpois, parName="lambda" )

# Posterior descriptives:
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples_genpois[,"lambda"] , main="lambda" , xlab=bquote(lambda) )

# Re-plot with different annotations:
plotPost( codaSamples_genpois[,"lambda"] , main="lambda" , xlab=bquote(lambda) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )

# Convergence diagnostics:
diagMCMC(codaObject=codaSamples_genpois, parName="omega" )

# Posterior descriptives:
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples_genpois[,"omega"] , main="omega" , xlab=bquote(omega) )

# Re-plot with different annotations:
plotPost( codaSamples_genpois[,"omega"] , main="omega" , xlab=bquote(omega) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )

# ---- 3. Data and model visualization ---------------------------- 
df_pois = with(data.frame(y_long), data.frame(x = y_long, y = dpois(y_long, m1_lambda)))
df_negbin = with(data.frame(y_long), data.frame(x = y_long, y = dnbinom(y_long, 
                                                                        m2_theta, 
                                                                        prob = m2_theta/(m2_theta+m2_lambda))))
df_genpois = with(data.frame(y_long), data.frame(x = y_long, y = dgpois(y_long, 
                                                                        m3_lambda, 
                                                                        m3_omega)))

ggplot(data.frame(y_long)) +
  geom_histogram(aes(x = y_long, y = ..density..),
                 binwidth = 1, fill = "grey", color = "black") +
  geom_line(data = df_pois, aes(x = x, y = y, color = 'red'), 
            linetype = 'solid', show.legend = T, size = 1) +
  geom_line(data = df_negbin, aes(x = x, y = y, color = "green"), 
            linetype = 'longdash', show.legend = T, size = 1) +
  geom_line(data = df_genpois, aes(x = x, y = y, color = "blue"), 
            linetype = "dotdash", show.legend = T, size = 1) +
  scale_color_discrete(name = "Models", 
                       labels = c(paste("M1: Poisson (", round(m1_lambda, 2), ")"), 
                                  paste("M2: Negative binomial (", round(m2_lambda, 2), ", ", round(m2_theta, 2), ")"), 
                                  paste("M3: Generalized Poisson (", round(m3_lambda, 2), ", ", round(m3_omega, 2), ")"))) + 
  labs(#title = titleM3,
       x = "# of claims",
       y = 'Share of claims / probability density') 
