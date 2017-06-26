################################################################################
##
## Build State Space Model by Rstan.
## Simulate ROI data & estimate time-variant regression coefficient
##
## written by Y.Nakahashi 
## 2017/6/25
##
################################################################################

################################################################################
## Environmental Settings
################################################################################
## set working directory
work_dir <- "/Users/nakahashi/Desktop/GitTest/MarketingMixModeling" # Private PC
setwd(work_dir)

## load libraries
library(dplyr)
library(tidyr)
library(tibble)

## load rstan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


################################################################################
## Data Simulation
################################################################################
####### simulation parameters settings
set.seed(123)
simulate_y <- function(pars) {
   n         <- pars[1] # num of observation
   mu        <- pars[2] # intercept
   var_e     <- pars[3] # residual variance
   beta_01   <- pars[4] # regression coefficient of X1 to be esitmated
   lambda_01 <- pars[5] # decay rate of Ad-Stock effect of X1
   beta_02   <- pars[6] # regression coefficient of X2 to be esitmated
   lambda_02 <- pars[7] # decay rate of Ad-Stock effect of X2
   
   X_01_raw <- rgamma(n, 3) * ifelse(runif(n) > 0.7, 0, 1)
   X_01_fil <- stats::filter(X_01_raw, lambda_01, "recursive")
   
   X_02_raw <- rgamma(n, 2) * ifelse(runif(n) > 0.8, 0, 1)
   X_02_fil <- stats::filter(X_02_raw, lambda_02, "recursive")
   
   error <- rnorm(n, 0, sqrt(var_e))
   
   y     <- mu + beta_01 * X_01_fil + beta_02 * X_02_fil + error
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

init_par <- array(c(100, 5, 2, 0.5, 0.6, 0.8, 0.5))
dat_Ana  <- na.omit(simulate_y(init_par))

################################################################################
## Run stan
################################################################################
## model
dat_Ord <- list(N    = nrow(dat_Ana),
                Y    = dat_Ana$Y,
                # K = 1,
                X_01 = dat_Ana$X_01,
                X_02 = dat_Ana$X_02)

## fitting 
fit_01 <- stan(file = './HierarchicalSSM_TimeVariantRegCoef_Sim.stan', 
               data = dat_Ord, 
               iter = 10000,
               chains = 4,
               seed = 1234,
               warmup = 5000)


################################################################################
## extract results
################################################################################
## Sample
res_01 <- rstan::extract(fit_01)

## parameters
ests <- summary(fit_01)$summary
beta_rows   <- rownames(ests)[grep("beta_*", rownames(ests))]
lambda_rows <- rownames(ests)[grep("lambda_*", rownames(ests))]

## Plot
pars <- c("mu", "var_error", "beta_X_01", "lambda_X_01", "beta_X_02",
          "lambda_X_02")
print(fit_01)
stan_trace(fit_01, pars = pars)
stan_hist(fit_01, pars = pars)
stan_dens(fit_01, separate_chains = T, pars = pars)
stan_ac(fit_01, separate_chains = T, pars = pars)

## Accuracy
data_stats <- c(mean(dat_Ana$Y), var(dat_Ana$True_Error), init_par[-c(1:3)])
Comp_Parameters <- data.frame(
   "True_Parameter" = c(init_par[-1]),
   "Data_Parameter" = data_stats,
   "Estimate"  = c(ests[1:length(pars), 1])) %>%
   rownames_to_column(var = "Parameters")



################################################################################
## Create output
################################################################################
## Regression coefficient
beta_par <- ests %>% 
   data.frame %>% 
   rownames_to_column(var = "Parameters") %>% 
   filter(Parameters %in% beta_rows)

## Decay rate
lambda_par <- ests %>% 
   data.frame %>% 
   rownames_to_column(var = "Parameters") %>% 
   filter(Parameters %in% lambda_rows)