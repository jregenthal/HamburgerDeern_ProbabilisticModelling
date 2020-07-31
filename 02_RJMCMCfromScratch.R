
# ---- Probabilistic Modelling ---------------------------------------------

# Author: Hamburger Deern

# Goal: Reversible Jump MCMC
## M1: Poisson
## M2: Negative binomial
## M3: Generalized poisson


# ---- 0. Preliminaries ---------------------------------------------------
#install.packages("LaplacesDemon")
#install.packages("runjags")
#install.packages("coda")
#install.packages("rjags")


library(LaplacesDemon)
library(runjags)
library(coda)
library(rjags)
library(dplyr)



# ---- 1. Data --------------------------------------------------------------------

#load data:
simulatedData <- read.csv("simulatedData.csv", row.names = 1)

y_long <- as.integer(simulatedData$sim_negbin)


# ---- 2. Functions ---------------------------------------------------------------

##joined probabilities of the models:
JP_m1 <- function(la) {
  P_y <- sum(dpois(y_long,la))
  P_lambda <- dgamma(la,1,1)
  P_m1 <- 1/3
  z <- P_y*P_lambda*P_m1
  return(z)
}

JP_m2 <- function(la, the) {
  p <- the/(la+the)
  P_y <- sum(dnbinom(y_long,the,p))
  P_theta <- 1/2*la*the^(-2)*(1+la/the)^(-3/2)
  P_lambda <- dgamma(la,1,1)
  P_m2 <- 1/3
  z <- P_y*P_theta*P_lambda*P_m2
  return(z)
}

JP_m3 <- function(la, ome) {
  P_y <- sum(dgpois(y_long,la,ome))
  the <- la*(1-ome)^2/(ome*(2-ome))
  P_omega <- 1/2*la*the^(-2)*(1+la/the)^(-3/2)
  P_lambda <- dgamma(la,1,1)
  P_m3 <- 1/3
  z <- P_y*P_omega*P_lambda*P_m3
  return(z)
}

#MCMC of the models:
MCMC_run <- function(model){
  #function for MCMC run
  #input: model1,2,3
  #output: mean of the computed parameters for the corresponding model
  if (newModel==1){
    ####MCMC of poisson model
    #Define model:
    modelString <- "
    model{
    for (i in 1:N) {
    y[i] ~ dpois(lambda)
    }
    lambda ~ dgamma(1,1)
    }
    "
    writeLines(modelString, con="claim_model.txt" )
    
    ##compile model:
    dataList <- list("y" = y_long, "N" = length(y_long))
    
    jagsModel <- jags.model(file="claim_model.txt",
                            data = dataList,
                            n.chains=1,
                            n.adapt=n.adapt)
    ##burn in:
    update(jagsModel, n.iter=n.burn)
    ##sample parameter:
    codaSamples = coda.samples(jagsModel, variable.names=c("lambda"), n.iter=n.iter)
    summary_CS <- summary(codaSamples)
    z <- as.numeric(summary_CS[[1]][1])
    return(z)
  }else if (newModel==2){
    
    ####MCMC of negative binomial model
    #Define model:
    
    modelString <- "
    model{
    for (i in 1:N) {
    y[i] ~ dnegbin(p[i],theta)
    p[i] <- theta/(theta+lambda)
    }
    theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
    omega ~ dbeta(omegaInit/(1-omegaInit),1)
    lambda ~ dgamma(1,1)
    }
    "
    writeLines(modelString, con="claim_model.txt" )
    
    ##compile model:
    omegaInit <- omega[length(omega)]
    dataList <- list("y" = y_long, "N" = length(y_long), "omegaInit" = omegaInit )
    
    jagsModel <- jags.model(file="claim_model.txt",
                            data = dataList,
                            n.chains=1,
                            n.adapt=n.adapt)
    ##burn in:
    update(jagsModel, n.iter=n.burn)
    
    ##sample parameter:
    codaSamples = coda.samples(jagsModel, variable.names=c("lambda", "theta"), n.iter=n.iter)
    summary_CS <- summary(codaSamples)
    x <- as.numeric(summary_CS[[1]][1])
    y <- as.numeric(summary_CS [[1]][2])
    z <- c(x,y)
    return(z)
  }else {
    ####MCMC of generalized poisson model
    #Define model:
    modelString <- "
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
    omega ~ dbeta(omegaInit/(1-omegaInit),1)
    lambda ~ dgamma(1,1)
    }
    "
    writeLines(modelString, con="claim_model.txt" )
    
    ##compile model:
    omegaInit <- omega[length(omega)]
    dataList <- list("y" = y_long, "N" = length(y_long), "omegaInit" = omegaInit )
    #initsList <- list(b=bDash_q_omega_m1)
    
    jagsModel <- jags.model(file="claim_model.txt",
                            data = dataList,
                            n.chains=1,
                            n.adapt=n.adapt)
    ##burn in:
    update(jagsModel, n.iter=n.burn)
    
    ##sample parameter:
    codaSamples = coda.samples(jagsModel, variable.names=c("lambda", "omega"), n.iter=n.burn)
    summary_CS <- summary(codaSamples)
    x <- as.numeric(summary_CS[[1]][1])
    y<- as.numeric(summary_CS[[1]][2])
    z <- c(x,y)
    return(z)
  }
}



# ---- 3. Pilot run ---------------------------------------------------------------
###Define parameters
n <- length(y_long) #number of rows
n.adapt = 100
n.burn = 100
n.iter = 1000 


####MCMC of negative binomial model
#Define model:
modelString <- "
model{
    for (i in 1:N) {
      y[i] ~ dnegbin(p[i],theta)
      p[i] <- theta/(theta+lambda)
    }
    theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
    omega ~ dbeta(1,1)
    lambda ~ dgamma(1,1)
}
"
writeLines(modelString, con="claim_model.txt" )

##compile model:
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel <- jags.model(file="claim_model.txt",
                        data = dataList,
                        n.chains=1,
                        n.adapt=n.adapt)
##burn in:
update(jagsModel, n.iter=n.burn)

##sample parameter:
codaSamples_NB = coda.samples(jagsModel, variable.names=c( "theta"), n.iter=n.iter)

####MCMC of generalized poisson model
#Define model:
modelString <- "
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
      theta <- lambda*(1-omega)*(1-omega)/(omega*(2-omega))
      omega ~ dbeta(1,1)
      lambda ~ dgamma(1,1)
    }
    "
writeLines(modelString, con="claim_model.txt" )

##compile model:
dataList <- list("y" = y_long, "N" = length(y_long))

jagsModel <- jags.model(file="claim_model.txt",
                        data = dataList,
                        n.chains=1,
                        n.adapt=n.adapt)
##burn in:
update(jagsModel, n.iter=n.burn)

##sample parameter:
codaSamples_GP = coda.samples(jagsModel, variable.names=c( "omega"), n.iter=n.iter)


#calculate mean and variance for q(v|m1) and q(w|m1): 
summary_CS_GP <- summary(codaSamples_GP)
omega_CS_mean <- as.numeric(summary_CS_GP[[1]][1])
omega_CS_var <- as.numeric(summary_CS_GP[[1]][2])^2
#a=mean(w)*((mean(w)*(1-mean(w))/var(w))-1) & b=a*(1-mean(w))/mean(w)
omega_mean <- omega_CS_mean*((omega_CS_mean*(1 - omega_CS_mean)/omega_CS_var)-1)
omega_var <- omega_mean*(1-omega_CS_mean)/omega_CS_mean


summary_CS_log_NB <- summary(codaSamples_NB)
theta_mean <- as.numeric(summary_CS_log_NB[[1]][1])
theta_var <- as.numeric(summary_CS_log_NB[[1]][2])^2



#---- 4. Initialization ----------------------------------------------------------



###Define parameters
IterationRJMCMC <- 250
BurnIn <- 51 # 50+1, because first row in parameter data frame will be the initialized values and needs to be excluded

###initialize start parameters
lambda_0 <- mean(y_long)
theta_0 <- max(c(0.01,mean(y_long)^2/(var(y_long) - mean(y_long))))
omega_0 <- max(c(0.01,1-sqrt(mean(y_long)/var(y_long))))

##vectors with all values for each parameter
lambda <- lambda_0
theta <- theta_0
omega <- omega_0

#data frame containing all calculated parameters of all MCMC iterations:
Parameters <- data.frame(matrix(0, ncol = 4, nrow = IterationRJMCMC+1))
colnames(Parameters) <-  c("lambda","theta","omega","model")
Parameters[1,] <- c(lambda_0,theta_0,omega_0,NA)



#---- 5. RJMCMC ---------------------------------------------------------

#### MCMC between models:
##choose first model to start with:
beta <- runif(1)
if (beta<=0.33){
  newModel <- 1
} else if (0.33<beta && beta<=0.66) {
  newModel <- 2
} else{
  newModel <- 3
}


##Start with RJMCMC:

for(i in 1:IterationRJMCMC){
  ##1: MCMC inside models:
  if (newModel==1){
    #MCMC:->newlambda
    newlambda <- MCMC_run(newModel)
    Parameters[i+1,1] <- newlambda
    Parameters[i+1,2] <- 0
    Parameters[i+1,3] <- 0
    Parameters[i+1,4] <- newModel
    lambda <- c(lambda, newlambda)
  } else if (newModel==2){
    #MCMC:->newlambda
    #MCMC:->newtheta
    MCMC_result <- MCMC_run(newModel)
    newlambda <-  MCMC_result[1]
    newtheta <-  MCMC_result[2]
    Parameters[i+1,1] <- newlambda
    Parameters[i+1,2] <- newtheta
    Parameters[i+1,3] <- 0
    Parameters[i+1,4] <- newModel
    lambda <- c(lambda, newlambda)
    theta <- c(theta, newtheta)
  } else {
    #MCMC:->newlambda
    #MCMC:->newomega
    MCMC_result <- MCMC_run(newModel)
    newlambda <- MCMC_result[1]
    newomega <- MCMC_result[2]
    Parameters[i+1,1] <- newlambda
    Parameters[i+1,2] <- 0
    Parameters[i+1,3] <- newomega
    Parameters[i+1,4] <- newModel
    lambda <- c(lambda, newlambda)
    omega <- c(omega, newomega)
  }
  oldModel <- newModel
  
  
  ##2: Propose new model to jump to:
  beta <- runif(1)
  if (oldModel==1){
    if (beta<0.5){
      proposedModel <- 2
    } else{
      proposedModel <- 3
    }
  }else if (oldModel==2){
    if (beta<0.5){
      proposedModel <- 1
    } else{
      proposedModel <- 3
    }
  }else{
    if (beta<0.5){
      proposedModel <- 1
    } else{
      proposedModel <- 2
    }
  }

  
  
  ##3: Calculate additional parameter:
  if (oldModel==1){
    if(proposedModel==2){
      thetaJump <- rlnorm(1, theta_mean,  theta_var)
    } else if (proposedModel==3){
      omegaJump <- rbeta(1, omega_mean,omega_var)
    }
  }else if (oldModel==2 && proposedModel==3){
    omegaJump <- 1-1/sqrt((1+lambda[length(lambda)]*1/theta[length(theta)]))
  }else if (oldModel==3 && proposedModel==2){
    thetaJump <- lambda[length(lambda)]*(1-omega[length(omega)])^2/(omega[length(omega)]*(2-omega[length(omega)]))
  }
  
  
  ##4: Calculate acceptance probability of model:
  if (oldModel==1){
    if (proposedModel==2){
      q_theta <- dlnorm(thetaJump,theta_mean, theta_var)
      delta <-  JP_m2(lambda[length(lambda)], thetaJump)/(JP_m1(lambda[length(lambda)]))*q_theta
      alpha <- min(1,delta)
    }else{
      #proposedModel==3
      q_omega <- dbeta(omegaJump,omega_mean,omega_var)
      delta <- JP_m3(lambda[length(lambda)], omegaJump)/(JP_m1(lambda[length(lambda)]))*q_omega
      alpha <- min(1,delta)
    }
  }else if (oldModel==2){
    if (proposedModel==1){
      q_theta <- dlnorm(theta[length(theta)],theta_mean, theta_var)
      delta <-  JP_m1(lambda[length(lambda)])/(JP_m2(lambda[length(lambda)], theta[length(theta)])*q_theta)
      alpha <- min(1,delta)
    }else{
      #proposedModel==3
      delta <- JP_m3(lambda[length(lambda)], omegaJump)/JP_m2(lambda[length(lambda)],theta[length(theta)])
      alpha <- min(1,delta)
    }
  }else {
    #oldModel==3
    if (proposedModel==1){
      q_omega <- dbeta(omega[length(omega)],omega_mean,omega_var)
      delta <- JP_m1(lambda[length(lambda)])/(JP_m3(lambda[length(lambda)], omega[length(omega)])*q_omega)
      alpha <- min(1,delta)
    }else{
      #proposedModel==2
      delta <- JP_m2(lambda[length(lambda)],thetaJump)/(JP_m3(lambda[length(lambda)], omega[length(omega)]))
      alpha <- min(1,delta)}
  }
  
  
  ##Decide if to jump to proposed model & save new parameters:
  u <- runif(1)
  if (alpha > u){
    newModel <- proposedModel
    if (oldModel==3 && newModel==2){
      theta <- c(theta,thetaJump)
    }else if (oldModel==2 && newModel==3){
      omega <- c(omega, omegaJump)
    }else if (oldModel==1 && newModel==2){
      theta <- c(theta, thetaJump)
    }else if (oldModel==1 && newModel==3){
      omega <- c(omega,omegaJump)
    }
  }
  
  
}


#---- 6.Model Diagnostics -------------------------------------------------------
#amount of jumps we have done inside each model:
Jumps_m1 <- as.numeric(count(filter(Parameters[BurnIn+1:length(Parameters),], model==1)))
Jumps_m2 <- as.numeric(count(filter(Parameters[BurnIn+1:length(Parameters),], model==2)))
Jumps_m3 <- as.numeric(count(filter(Parameters[BurnIn+1:length(Parameters),], model==3)))
#posterior model probability:
posterioModelProb_m1 <- 1/(IterationRJMCMC-BurnIn)*Jumps_m1
posterioModelProb_m2 <- 1/(IterationRJMCMC-BurnIn)*Jumps_m2
posterioModelProb_m3 <- 1/(IterationRJMCMC-BurnIn)*Jumps_m3
#posterior model odds:
PO_m1_m2 <- posterioModelProb_m1/posterioModelProb_m2
PO_m1_m3 <- posterioModelProb_m1/posterioModelProb_m3
PO_m2_m3 <- posterioModelProb_m2/posterioModelProb_m3
PO_m2_m1 <- posterioModelProb_m2/posterioModelProb_m1
PO_m3_m1 <- posterioModelProb_m3/posterioModelProb_m1
PO_m3_m2 <- posterioModelProb_m3/posterioModelProb_m2
#log bayes factor: =log(PO) since f(m)=1/2 for all m
BF_12 <- log(posterioModelProb_m1/posterioModelProb_m2)
BF_13 <- log(posterioModelProb_m1/posterioModelProb_m3)
BF_23 <- log(posterioModelProb_m2/posterioModelProb_m3)
BF_21 <- log(posterioModelProb_m2/posterioModelProb_m1)
BF_31 <- log(posterioModelProb_m3/posterioModelProb_m1)
BF_32 <- log(posterioModelProb_m3/posterioModelProb_m2)


