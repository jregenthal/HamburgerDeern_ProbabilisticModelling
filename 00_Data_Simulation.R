# ---- Probilistic Modelling --------------------------------------
#
# Author: Hamburger Deern
# 
# Goal: Simulate three datasets from distributions
# - Poisson
# - Negative Binomial
# - Generalized Poisson


# ---- 0. Preliminaries -------------------------
library(ggplot2)
library(HMMpa)

set.seed(1)

lambda = rgamma(1, 1, 1)
omega = rbeta(1, 1, 1)
theta = lambda*(1-omega)*(1-omega)/(omega*(2-omega))


# ---- 1. Simulation -------------------------
datapoints = 100000 

## M1: Poisson
sim_pois = rpois(datapoints, lambda)

ggplot(data = data.frame(sim_pois), aes(x = sim_pois)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black")

## M2: Negative binomial
sim_negbin = rnbinom(datapoints, theta, 
                     prob = theta/(lambda + theta))
ggplot(data = data.frame(sim_negbin), aes(x = sim_negbin)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black")


## M3: Generalized Poisson
# documentation HMMpa: https://www.rdocumentation.org/packages/HMMpa/versions/1.0.1/topics/dgenpois
# dgpois(c(1,2,3), lambda, omega)
# dgenpois(c(1,2,3), (1-omega)*lambda, omega)
sim_genpois = rgenpois(datapoints, lambda1 = (1-omega)*lambda, lambda2 = omega) 
ggplot(data = data.frame(sim_genpois), aes(x = sim_genpois)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black")


simulatedData <- cbind(sim_pois, sim_negbin, sim_genpois)
write.csv(simulatedData, "00_Data_Simulation.csv")
