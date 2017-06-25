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

## load rstan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


################################################################################
## Data Simulation
## Assume following model:
##   Model 1 : log(Y_ti) ~ state_ti + beta_ti * X_ti + error_ti
################################################################################
####### simulation parameters settings
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












## set seed
set.seed(123)

## number of data
num_week    <- 12
num_year    <- 4
num_region  <- 5
data_length <- num_week * num_year

## distributions of
## state_t0
mu_state_t0  <- 3
var_state_t0 <- 0.1

## state_t
var_state <- 0.01

## white noise
var_error <- 0.000025

## time-variant regression coefficients of X
beta_X_0   <- sin(seq(1, 3, length=data_length))/100
var_beta_X <- 0.000001



##### simulate data
## X
X_0 <- ceiling(cos(5:17)*10+10)
X   <- matrix(numeric(data_length * num_region), nrow = data_length)
scl_TV <- 5

for (i in 1:num_region) {
   for (j in 1:num_year) {
      for(k in 1:num_week) {
         if (j == 1) {
            X[(j-1)*12+k, i] <- X_0[k] + rgamma(1, scl_TV)
         } else {
            X[(j-1)*12+k, i] <- X[(j-2)*12+k, i] + 
               rgamma(1, scl_TV)
         }
      }
   }
}
X_DF <- as.data.frame(X) %>% 
   gather() %>% 
   mutate("Area" = gsub("V", "", key)) %>%
   select(Area, value)


## state_t
state      <- matrix(numeric(length = data_length * num_region), 
                     nrow = data_length)
state[1, ] <- rnorm(num_region, mean = mu_state_t0, sd = sqrt(var_state_t0))

for (i in 1:num_region) {
   for(j in 2:data_length) {
      state[j, i] <- rnorm(1, mean = state[j-1, i], sd = sqrt(var_state))
   }
}
state_DF <- as.data.frame(state) %>% 
   gather() %>% 
   mutate("Area" = gsub("V", "", key)) %>%
   select(Area, value)


## beta of X
beta_X <- matrix(numeric(data_length * num_region), nrow = data_length)
beta_X[1, ] <- 
for (i in 1:num_region) {
   for (j in 2:data_length) {
      beta_X[j, i] <- beta_X_0[j] + rnorm(1, 0, sd = sqrt(var_TV))
   }
}
beta_X_DF <- as.data.frame(beta_X) %>% 
   gather() %>% 
   mutate("Area" = gsub("V", "", key)) %>%
   select(Area, value)


## white noise
error <- matrix(rnorm(data_length * num_region, 0, sqrt(var_error)), 
                nrow=data_length)
error_DF <- as.data.frame(error) %>% 
   gather() %>% 
   mutate("Area" = gsub("V", "", key)) %>%
   select(Area, value)


## create data
log_y <- state_DF$value + beta_X_DF$value * X_DF$value + error_DF$value
Y     <- ceiling(exp(log_y))
log_y_real <- log(Y)

dat_Ana <- data.frame(
   "OrderNum" = Y,
   "LogOrder" = log_y_real,
   "Area"     = X_DF$Area,
   "X"        = X_DF$value
)

# plot(dat_Ana$OrderNum, type="l")
# plot(dat_Ana$LogOrder, type="l")
# plot(dat_Ana$X, type="l")

################################################################################
## Run stan
################################################################################
## model
dat_Ord <- list(N = nrow(dat_Ana)/length(unique(dat_Ana$Area)),
                Y = dat_Ana$LogOrder,
                K = length(unique(dat_Ana$Area)),
                X = dat_Ana$X)

## fitting 
fit_01 <- stan(file = './StanModel/HierarchicalSSM_TimeVariantRegCoef_Sim.stan', 
               data = dat_Ord, 
               iter = 10000,
               chains = 4,
               seed = 1234)


################################################################################
## extract results
################################################################################
## sample
res_01 <- rstan::extract(fit_01)

## parameters
ests <- summary(fit_01)$summary
state_rows <- rownames(ests)[grep("state_t*", rownames(ests))]
beta_rows  <- rownames(ests)[grep("beta_*", rownames(ests))]

## state
state_par <- ests %>% data.frame %>% select(mean) %>% 
   mutate("Par" = rownames(ests)) %>% 
   filter(Par %in% state_rows)
est_state <- state_par$mean[2:(length(state_par$mean)-1)]
plot(est_state, type="l")


## regression coefficient
beta_par <- ests %>% data.frame %>% select(mean) %>% 
   mutate("Par" = rownames(ests)) %>% 
   filter(Par %in% beta_rows)
est_beta <- beta_par$mean[2:(length(beta_par$mean)-1)]
plot(est_beta, type="l")


## compare true & estimate parameters
plot(cbind(state_DF$value, est_state))
plot(state_DF$value, type="b")
lines(est_state, col=2) 

plot(cbind(beta_X_DF$value, est_beta))
plot(beta_X_DF$value, type="b")
lines(est_beta, col=2) 

y_hat <- est_state + dat_Ana$X * est_beta
plot(dat_Ana$LogOrder, y_hat)
plot(dat_Ana$LogOrder, type="b")
lines(y_hat, col=2)

## sampling check
stan_trace(fit_01)
stan_hist(fit_01)
stan_dens(fit_01, separate_chains = T)
stan_ac(fit_01, separate_chains = T)
