################################################################################
##
## Estimate Ad-stock effect.
## written by Y.Nakahashi 
## 2017/6/2
##
################################################################################

################################################################################
## Environmental Settings
################################################################################
## set working directory
work_dir <- "/Users/nakahashi/Desktop/GitTest/MarketingMixModeling"
setwd(work_dir)

## libraries
library(dplyr)
library(tidyr)

################################################################################
## Palda model
################################################################################
simulate_y <- function(pars) {
   n         <- pars[1] # num of observation
   mu        <- pars[2] # intercept
   beta_01   <- pars[3] # regression coefficient of X1 to be esitmated
   beta_02   <- pars[4] # regression coefficient of X2 to be esitmated
   lambda_01 <- pars[5] # decay rate of Ad-Stock effect of X1
   lambda_02 <- pars[6] # decay rate of Ad-Stock effect of X2
   var_e     <- pars[7] # residual variance
   
   X_01_raw <- rnorm(n, 100, 2) * ifelse(runif(n) > 0.7, 0, 1)
   X_01_fil <- stats::filter(X_01_raw, lambda_01, "recursive")

   X_02_raw <- rnorm(n, 100, 2) * ifelse(runif(n) > 0.2, 0, 1)
   X_02_fil <- stats::filter(X_02_raw, lambda_02, "recursive")
   
   error <- rnorm(n, 0, sqrt(var_e))

   y     <- mu + beta_01 * X_01_fil + beta_02 * X_02_fil + error
   dat <- data.frame(
      "Y" = y,
      "X_01" = X_01_raw,
      "X_02" = X_02_raw,
      "Y_lag" = dplyr::lag(y, 1))
   return(na.omit(dat))
}

## check
init_par <- array(c(100, 5, 0.5, 0.3, 0.9, 0.2, 2))
simulate_y(init_par)


## simulation
set.seed(123)
init_par      <- array(c(100, 5, 0.5, 0.3, 0.9, 0.2, 5))
k             <- 500
results_palda <- matrix(0, k, 4)
for (i in 1:k) {
   dat <- simulate_y(init_par)
   res <- lm(Y ~ X_01 + X_02 + Y_lag, data = dat)
   results_palda[i, ] <- coef(res)
}

summary(results_palda)
MASS::truehist(results_palda[, 1])
MASS::truehist(results_palda[, 2])
MASS::truehist(results_palda[, 3])
MASS::truehist(results_palda[, 4])


################################################################################
## Optim
################################################################################
return_AIC <- function(param) {
   a <- param[1]
   b <- param[2]
   dat$X_01_fil <- stats::filter(dat$X_01, a, "recursive")
   dat$X_02_fil <- stats::filter(dat$X_02, b, "recursive")
   AIC(lm(Y ~ X_01_fil + X_02_fil, dat))
}

## check
init_par  <- array(c(100, 5, 0.5, 0.3, 0.9, 0.2, 2))
dat       <- simulate_y(init_par)
decay_par <- array(c(0.5, 0.5))
optim(par = optim(par = decay_par, fn = return_AIC)$par, 
      fn = return_AIC)$par


set.seed(456)
init_par      <- array(c(100, 5, 0.5, 0.3, 0.9, 0.2, 5))
decay_par     <- array(c(0.5, 0.5))
k             <- 500
results_optim <- matrix(0, k, 5)
for (i in 1:k) {
   dat <- simulate_y(init_par)
   res <- optim(par = optim(par = decay_par, fn = return_AIC)$par, 
                fn = return_AIC)$par
   
   dat$X_01_fil <- stats::filter(dat$X_01, res[1], "recursive")
   dat$X_02_fil <- stats::filter(dat$X_02, res[2], "recursive")

   lm_out             <- lm(Y ~ X_01_fil + X_02_fil, data = dat)
   results_optim[i, ] <- c(res, coef(lm_out))
   if (i %% 50 == 0) cat("Round: ", i, "\n")
}

summary(results_optim)
MASS::truehist(results_optim[, 1])
MASS::truehist(results_optim[, 2])
MASS::truehist(results_optim[, 3])
MASS::truehist(results_optim[, 4])



################################################################################
## GA
################################################################################
# install.packages("GA")
library(GA)

## Data simulation function
simulate_y <- function(pars) {
   ## Simulation parameters
   n         <- pars[1] # num of observation
   mu        <- pars[2] # intercept
   var_e     <- pars[3] # residual variance
   beta_01   <- pars[4] # regression coefficient of X1 to be esitmated
   lambda_01 <- pars[5] # decay rate of Ad-Stock effect of X1
   beta_02   <- pars[6] # regression coefficient of X2 to be esitmated
   lambda_02 <- pars[7] # decay rate of Ad-Stock effect of X2
   
   ## Create true Ad-Stock variables
   X_01_raw <- rgamma(n, 3) * ifelse(runif(n) > 0.7, 0, 1)
   X_01_fil <- stats::filter(X_01_raw, lambda_01, "recursive")
   
   X_02_raw <- rgamma(n, 2) * ifelse(runif(n) > 0.8, 0, 1)
   X_02_fil <- stats::filter(X_02_raw, lambda_02, "recursive")
   
   ## Create residuals
   error <- rnorm(n, 0, sqrt(var_e))
   
   ## Create observations   
   y     <- mu + beta_01 * X_01_fil + beta_02 * X_02_fil + error
   
   ## Return dataset
   dat <- data.frame(
      "Y"          = y,
      "X_01"       = X_01_raw,
      "X_02"       = X_02_raw,
      "X_01_Fil"   = X_01_fil,
      "X_02_Fil"   = X_02_fil,
      "Y_lag"      = dplyr::lag(y, 1),
      "True_Error" = error)
   return(dat)
}

## Data simulation
set.seed(123)
init_par <- array(c(100, 5, 2, 0.5, 0.6, 0.8, 0.5))
dat_Ana  <- na.omit(simulate_y(init_par))

## Define objective function
OLS <- function(dat, pars) {
   X_01_Fil_E <- stats::filter(dat$X_01, pars[3], "recursive")
   X_02_Fil_E <- stats::filter(dat$X_02, pars[5], "recursive")
   Y_hat <- pars[1] + pars[2] * X_01_Fil_E + pars[4] * X_02_Fil_E
   SSE   <- t(dat$Y - Y_hat) %*% (dat$Y - Y_hat)
   return(SSE)
}


## Fit via GA
res_GA <- ga(type = 'real-valued', 
             min = c(0, 0, 0, 0, 0), 
             max = c(20, 2, 1, 2, 1),
             popSize = 500, maxiter = 500, 
             names = c('Intercept', 'b_01', 'l_01', 'b_02', 'l_02'),
             keepBest = T, 
             fitness = function(b) -OLS(dat_Ana, b))

## Result
summary(res_GA)$solution
init_par[c(-1, -3)]

-summary(res_GA)$fitness / (nrow(dat_Ana) - 6) ### Residual variance
init_par[3]



