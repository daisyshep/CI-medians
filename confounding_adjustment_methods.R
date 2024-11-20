############################################################################################################################################
## CONFOUNDING-ADJUSTMENT METHODS FOR THE DIFFERENCE IN MEDIANS
## Creator: D. A. Shepherd (UoM, MCRI)
## Collaborators: B. Baer (UoR), M. Moreno-Betancur (UoM, MCRI)
## Date created: August 2023
## Last updated: November 2024
############################################################################################################################################

## NOTES:
# Functions supporting the manuscript:
#       'Confounding-adjustment method for the difference in medians'
#       Published and available at: https://doi.org/10.1186/s12874-023-02100-6.
# Functions to estimate the causal difference in medians
# Bootstrap functions provided to estimate the std. errors and 95% CIs

## Please refer to the README for further details

## CITATION: Please use the following citations when using this resource:
# Shepherd, D., Baer, B. & Moreno-Betancur, M. Confounding-adjustment methods for the causal difference in medians. BMC Med Res Methodol 23, 288 (2023). https://doi.org/10.1186/s12874-023-02100-6.
# Shepherd, D., Baer, B. & Moreno-Betancur, M. CI-medians. https://github.com/daisyshep/CI-medians.

############################################################################################################################################

## DATA COMPATIBILITY: 
## Data file compatible with these functions is of the following form:
##       - Continuous outcome (denoted as Y)
##       - Binary exposure (denoted as A; coded as 0/1)
##       - Confounders (five in total; denoted as C1,...,C5)

######################################################################

## Loading required libraries (may require installation prior to use)
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
# G-COMPUTATION
# Two variations in implementation:
# 1. G-comp (MC)
# 2. G-comp (approx)
######################################################################

# G-COMPUTATION (MC)
gcomp_mc <- function(data, ind, out, exp, form, reps){

  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for gcompboot function
  #       - out = outcome variable
  #       - exp = exposure variable
  #       - form = formula for outcome model
  #       - reps = number of reps for sampling per observation (R, integer)

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
      
      ## Fitting the linear model
      fit <- lm(as.formula(formula), data=onesample)
      
      ## Obtain individual predictions under fitted model
      onesample$predicted_mY <- predict(fit, onesample)
      
      # Sampling to build exp and unexp distributions
      # Looping over number of obs
      obs_exp <- c()
      obs_unexp <- c()
      
      for (i in 1:nrow(dat)){
        obs_exp <- c(obs_exp,
                     rlnorm(reps, 
                            meanlog=onesample[which(onesample$index==1),]$predicted_mY,
                            sdlog=sigma(fit)))
        obs_unexp <- c(obs_unexp,
                       rlnorm(reps, 
                              meanlog=onesample[which(onesample$index==0),]$predicted_mY,
                              sdlog=sigma(fit)))
      } 
      
      ## Calculating average causal effect (ACE)
      med_exp <- median(obs_exp)
      med_unexp <- median(obs_unexp)
      ACE_diff_med <- med_exp - med_unexp
      return(ACE_diff_med)
    }
  }
  
  ## Not continuous outcome
  else{
    print("Error: Outcome must be continuous")
  } 
}


#######################################

## G-COMPUTATION (APPROX)
gcomp_approx <- function(data, ind, out, exp, formula, y_support){
  
  # ARGUMENTS:
  #       - data = data file as specified above
  #       - ind = argument required for gcompboot function
  #       - out = outcome variable
  #       - exp = exposure variable
  #       - form = formula for outcome model
  #       - y_support = vector of candidate y values (y*) to search over for the median
  #                     (Note: assume vector is equally spaced across range)
  
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
      
      ## Fitting the linear model
      fit <- lm(as.formula(formula), data=onesample)
      
      ## Obtain individual predictions under fitted model
      onesample$predicted_mY <- predict(fit, onesample)
      
      ## Looping over y vector and for exposed and unexposed
      g_y_exp <- c()
      g_y_unexp <- c()
      
      for (i in 1:length(y_support)){
        ## Estimate pdf for y|A=1,C
        p_y <- dlnorm(x=y_support[i],
                      meanlog=onesample[which(onesample$index==1),]$predicted_mY,
                      sdlog=sigma(fit))
        ## Take the mean of these values
        g_y_exp[i] <- mean(p_y)
        
        ## Estimate pdf for y|A=0,C
        p_y <- dlnorm(x=y_support[i],
                      meanlog=onesample[which(onesample$index==0),]$predicted_mY,
                      sdlog=sigma(fit))
        ## Take the mean of these values
        g_y_unexp[i] <- mean(p_y)
      }
      
      ## STEP 5: Obtain the estimate of the median
      increment <- y_support[2]-y_support[1]
      cdf_exp <- data.frame(y, g_y_cum=cumsum(increment*g_y_exp))
      cdf_unexp <- data.frame(y, g_y_cum=cumsum(increment*g_y_unexp))
      
      med_exp <- cdf_exp$y[which(cdf_exp$g_y_cum >= 0.5)[1]]
      med_unexp <- cdf_unexp$y[which(cdf_unexp$g_y_cum >= 0.5)[1]]
      
      ACE_diff_med <- med_exp - med_unexp
      return(ACE_diff_med) 
    }
  }
  
  ## Not continuous outcome
  else{
    print("Error: Outcome must be continuous")
  } 
}

#######################################
gcompboot <- function(runboot_dat, R, stat, out, exp, form, ncpu, ...){
  
  # ARGUMENTS:
  #       - runboot_dat = data file as specified above
  #       - R = number of bootstrap repetitions (calculating se and CIs)
  #       - stat = argument required for gcompboot function
  #       - form = formula for outcome model
  #       - ncpu = no. of processes to be used in parallel (integer)
  #       - ... = additonal arguments to be passed to stat function
  #               (e.g., reps for gcomp_sample; y_support for gcomp_grid)
  
  ## Performing bootstrap
  bstrap <- boot(data=runboot_dat, statistic=stat, stype="i", R=R,
                 out=out, exp=exp, formula=form,
                 parallel="multicore", ncpus=ncpu, ...)
  
  ## Obtaining summary statistics
  est <- bstrap$t0
  se <- apply(bstrap$t, 2, sd, na.rm=T)
  bt <- boot.ci(bstrap, index=1, type=c("perc"))
  ci_low <- bt$percent[4]
  ci_up <- bt$percent[5]
  return(c(est=est, se=se, ci_low=ci_low, ci_up=ci_up))
}
######################################################################
