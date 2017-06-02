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

library(dplyr)

################################################################################
## Palda model
################################################################################
simulate_y <- function(pars) {
   n      <- pars[1] # num of observation
   mu     <- pars[2] # intercept
   beta   <- pars[3] # regression coefficient to be esitmated
   lambda <- pars[4] # decay rate of Ad-Stock effect
   var_e  <- pars[5] # residual variance
   
   X_raw <- rnorm(n, 100, 2) * ifelse(runif(n) > 0.7, 0, 1) # X
   X_fil <- stats::filter(X_raw, lambda, "recursive") # filter-ed X
   error <- rnorm(n, 0, sqrt(var_e)) # white noise
   y     <- mu + beta * X_fil + error # create y

   dat <- data.frame(
      "Y" = y,
      "X" = X_raw,
      "Y_lag" = dplyr::lag(y, 1))
   return(na.omit(dat))
}

## check
init_par <- array(c(100, 5, 0.5, 0.9, 2))
simulate_y(init_par)


## simulation
set.seed(123)
init_par <- array(c(100, 5, 0.5, 0.9, 2))
k        <- 500
results  <- matrix(0, k, 3)
for (i in 1:k) {
   dat <- simulate_y(init_par)
   res <- lm(Y ~ ., data = dat)
   results[i, 1:3] <- coef(res)
}

summary(results[, 2])
MASS::truehist(results[, 2])
summary(results[, 3])
MASS::truehist(results[, 3])



################################################################################
## Optim
################################################################################
simulate_y_02 <- function(pars) {
   X1 <- rnorm(100, 100, 2) * rep(c(0, rep(1, 7), 0, 1), 10)
   X2 <- filter(X1, pars[1], "recursive")
   Z1 <- rnorm(100, 5, 2) * rep(c(rep(0, 5), rep(1, 5)), 10)
   Z2 <- filter(Z1, pars[2], "recursive")
   dat <- data.frame(
      "Y"  = 50 + pars[3] * X2 + pars[4] * Z2 + rnorm(100, 0, pars[5]),
      "X1" = X1,
      "Z1" = Z1)
   return(dat)
}

return_AIC_02 <- function(param) {
   a <- param[1]
   b <- param[2]
   dat$TmpX <- filter(dat$X1, a, "recursive")
   dat$TmpZ <- filter(dat$Z1, b, "recursive")
   AIC(lm(Y ~ TmpX + TmpZ, dat))
}

init_par <- array(c(0.5, 0.5))
optim(par = optim(par = init_par, fn = return_AIC)$par, 
      fn = return_AIC)$par

pars      <- c(0.8, 0.7, 0.5, 4, 5.0)
n         <- 100
res_02_02 <- matrix(0, n, 5)
for (i in 1:n) {
   dat <- simulate_y(pars)
   res <- optim(par = optim(par = init_par, fn = return_AIC)$par, 
                fn = return_AIC)$par
   X2  <- filter(dat$X1, res[1], "recursive")
   Z2  <- filter(dat$Z1, res[2], "recursive")
   res_02_02[i, ] <- c(res, coef(lm(Y ~ X2 + Z2, data = dat)))
}
summary(res_02_02)


