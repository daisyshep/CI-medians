######################################################################
## CONFOUNDING-ADJUSTMENT METHODS - DIFFERENCE IN MEDIANS
## Creator: D. A. Shepherd (UoM, MCRI)
## Date edited: April 2022
######################################################################

# NOTES:
# Functions supporting the manuscript XXXX
# Functions to estimate the causal difference in medians
# Bootstrap funtions provided to estimate the std. errors and 95% CIs

# Data file compatible with these functions is of the following form:
#       - Continous outcome (denoted as Y)
#       - Binary exposure (denoted as A; coded as 0/1)
#       - Confounders (five in total; denoted as C1,...,C5)

######################################################################

## Loading required libraries
library(boot)
library(quantreg)
library(MASS)
library(Hmisc)


######################################################################
# QUANTILE REGRESSION:
# Outcome/quantile regression (no interactions; assumes constant CE)
######################################################################

outreg_med <- function(data, formula, R){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - formula = formula for outcome regression model
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  
  # Fitting the quantile regression model
  mod_rq <- rq(as.formula(formula), data=data, method="br")
  
  # Summary statistics
  est <- as.numeric(mod_rq$coef[2])
  
  # CI - bootstrap method
  bt <- boot.rq(cbind(1, data[,c("A","C1","C2","C3","C4","C5")]), data$Y, 
                tau=0.5, R=R, bsmethod="mcmb")
  se <- t(apply(bt$B, 2, sd))[2]
  ci_low <- t(apply(bt$B, 2, quantile, c(0.025)))[2]
  ci_up <- t(apply(bt$B, 2, quantile, c(0.975)))[2]
  
  return(c(est=est, se=se, ci_low=ci_low, ci_up=ci_up))
}
######################################################################



######################################################################
# WEIGHTED QUANTILE REGRESSION
# Quantile regression weighted via weights from propensity score model
######################################################################

weight_qr <- function(data, ind, form1, form2, R){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for weight_qr_boot function
  #       - form1 = formula for propensity score (PS) model
  #       - form2 = formula for outcome model conditional on exposure
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  
  dat <- data[ind,]
  
  # Step 1: Calculate the probability of exposure | covariates
  prob_mod <- glm(as.formula(form1), family=binomial(logit), data=dat)
  
  # Step 2: Calculate propensity scores
  ps <- predict(prob_mod, dat, type="response")

  # Step 3: Generate IP weights
  ip_weights <- ifelse(dat$A==1, 1/ps, 1/(1-ps))
  dat2 <- cbind(dat, ip_weights)
  
  # Step 4: Fit the quantile regression model
  mod_rq <- rq(as.formula(form2), tau=0.5, data=dat2, weights=ip_weights, method="br")
  
  # Getting summary stats
  est <- as.numeric(mod_rq$coef[2])
  return(est)
}
#######################################
weight_qr_boot <- function(runboot_dat, R, stat, form1, form2, ncpu){
  
  # ARGUMENTS:
  #       - runboot_dat = data file as specified above
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  #       - stat = argument required for weight_qr_boot function
  #       - form1 = formula for propensity score (PS) model
  #       - form2 = formula for outcome model conditional on exposure
  #       - ncpu = 
  
  ## Performing bootstrap
  bstrap <- boot(data=runboot_dat, statistic=stat, stype="i", R=R, form1=form1, form2=form2,
                 parallel="multicore", ncpus=ncpu)
  
  ## Obtaining summary statistics
  est <- bstrap$t0
  se <- apply(bstrap$t, 2, sd, na.rm=T)
  bt <- boot.ci(bstrap, index=1, type=c("perc"))
  ci_low <- bt$percent[4]
  ci_up <- bt$percent[5]
  return(c(est=est, se=se, ci_low=ci_low, ci_up=ci_up))
}
######################################################################



######################################################################
# IPW ESTIMATOR:
# Based on Zhang et al. (2012) paper (DOI: 10.1111/j.1541-0420.2011.01712.x)
######################################################################

ipw_estimator <- function(data, ind, form){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for weight_qr_boot function
  #       - form = formula for propensity score (PS) model
  
  dat <- data[ind,]
  
  # Step 1: Calculate weights
  ps_mod <- glm(as.formula(form), family=binomial(logit), data=dat)
  
  # Step 2: Calculate propensity scores
  ps <- predict(ps_mod, dat, type="response")

  # Step 3: Generate inverse probability weights + normalised
  ipw <- ifelse(dat$A==1, 1/ps, 1/(1-ps))
  ipw2 <- ifelse(dat$A==1, ipw/(sum(ipw[dat$A==1])), ipw/(sum(ipw[dat$A==0])))
  
  # Step 4: Function for uniroot
  f1 <- function(x, weights, outcome){
    (sum(weights*ifelse(outcome<=x, 1, 0)) - 0.5)
  }
  
  # Step 5: Solving for root
  rt_1 <- uniroot(f1, weights=ipw2[dat$A==1], outcome=dat$Y[dat$A==1], interval=c(-100,100))$root
  rt_0 <- uniroot(f1, weights=ipw2[dat$A==0], outcome=dat$Y[dat$A==0], interval=c(-100,100))$root
  return(rt_1-rt_0)
}
#######################################

ipw_boot <- function(runboot_dat, R, stat, form, ncpu){
  
  # ARGUMENTS:
  #       - runboot_dat = data file as specified above
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  #       - stat = argument required for bootstrap function (e.g., ipw_estimator)
  #       - form = formula for propensity score (PS) model
  #       - ncpu = no. of processes to be used in parallel (integer)
  
  ## Performing bootstrap
  bstrap <- boot(data=runboot_dat, statistic=stat, stype="i", R=R, form=form,
                 parallel="multicore", ncpus=ncpu)
  
  ## Obtaining summary statistics
  est <- bstrap$t0
  se <- apply(bstrap$t, 2, sd, na.rm=T)
  bt <- boot.ci(bstrap, index=1, type=c("perc"))
  ci_low <- bt$percent[4]
  ci_up <- bt$percent[5]
  return(c(est=est, se=se, ci_low=ci_low, ci_up=ci_up))
}
#######################################



######################################################################
# HEURISTIC G-COMPUTATION
# Two variations in implementation:
# 1. G-comp (med): Calculated as median of predictions
# 2. G-comp (mean): Calculated as mean of predictions
######################################################################

# G-COMP (MED)
gcompstat_med <- function(data, ind, out, exp, form){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for gcompboot function
  #       - out = outcome variable
  #       - exp = exposure variable
  #       - form = formula for outcome model

  dat <- data[ind,] 
  
  # Original dataset
  dat$index <- -1
  
  # Dataset untreated: Treatment set to 0
  dat0 <- dat
  dat0$index <- 0
  dat0[[exp]] <- 0
  dat0[[out]] <- NA
  
  # Dataset treated: Treatment set to 1
  dat1 <- dat
  dat1$index <- 1
  dat1[[exp]] <- 1
  dat1[[out]] <- NA
  
  # Combining into one sample
  onesample <- rbind(dat,dat0,dat1)
  
  ## Continuous Outcome
  if(is.factor(dat[[out]])==FALSE){
    if(is.numeric(dat[[out]])==TRUE){
      
      # Fitting the linear model
      fit <- rq(as.formula(form), data=onesample, method="br")
      
      ## Predicted values under fitted model
      onesample$predicted_medY <- predict(fit, onesample)
      
      ## Calculating average causal effect (ACE)
      ## Calculated as mean of medians
      EYA0 <- median(onesample[which(onesample$index==0),]$predicted_medY) 
      EYA1 <- median(onesample[which(onesample$index==1),]$predicted_medY) 
      ACE_diff_med <- EYA1 - EYA0
      return(ACE_diff_med)
    }
  }
  
  ## Not continuous outcome
  else{
    print("Error: Outcome must be continuous")
  } 
}

#######################################

# G-COMP (MEAN)

gcompstat_mean <- function(data, ind, out, exp, form){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for gcompboot function
  #       - out = outcome variable
  #       - exp = exposure variable
  #       - form = formula for outcome model
  
  dat <- data[ind,] 
  
  # Original dataset
  dat$index <- -1
  
  # Dataset untreated: Treatment set to 0
  dat0 <- dat
  dat0$index <- 0
  dat0[[exp]] <- 0
  dat0[[out]] <- NA
  
  # Dataset treated: Treatment set to 1
  dat1 <- dat
  dat1$index <- 1
  dat1[[exp]] <- 1
  dat1[[out]] <- NA
  
  # Combining into one sample
  onesample <- rbind(dat,dat0,dat1)
  
  ## Continuous Outcome
  if(is.factor(dat[[out]])==FALSE){
    if(is.numeric(dat[[out]])==TRUE){
      
      # Fitting the linear model
      fit <- rq(as.formula(form), data=onesample, method="br")
      
      ## Predicted values under fitted model
      onesample$predicted_medY <- predict(fit, onesample)
      
      ## Calculating average causal effect (ACE)
      ## Calculated as mean of medians
      EYA0 <- mean(onesample[which(onesample$index==0),]$predicted_medY) 
      EYA1 <- mean(onesample[which(onesample$index==1),]$predicted_medY) 
      ACE_diff_med <- EYA1 - EYA0
      return(ACE_diff_med)
    }
  }
  
  ## Not continuous outcome
  else{
    print("Error: Outcome must be continuous")
  } 
}

#######################################
gcompboot <- function(runboot_dat, R, stat, out, exp, form, ncpu){
  
  # ARGUMENTS:
  #       - runboot_dat = data file as specified above
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  #       - stat = argument required for gcompboot function
  #       - form = formula for outcome model
  #       - ncpu = no. of processes to be used in parallel (integer)
  
  ## Performing bootstrap
  bstrap <- boot(data=runboot_dat, statistic=stat, stype="i", R=R,
                 out=out, exp=exp, formula=form,
                 parallel="multicore", ncpus=ncpu)
  
  ## Obtaining summary statistics
  est <- bstrap$t0
  se <- apply(bstrap$t, 2, sd, na.rm=T)
  bt <- boot.ci(bstrap, index=1, type=c("perc"))
  ci_low <- bt$percent[4]
  ci_up <- bt$percent[5]
  return(c(est=est, se=se, ci_low=ci_low, ci_up=ci_up))
}
######################################################################