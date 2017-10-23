### library ##############################
library(GA)

### Parameters settings
n_year      <- 4
n_per_year  <- 12
n_per_month <- 4

n          <- n_year * n_per_year * n_per_month
mu         <- 4
var_e      <- 0.5

beta_01 <- 0.03
beta_02 <- 0.05
beta_03 <- 0.01
beta_04 <- 0.01
beta_05 <- 0.01

lambda_01 <- 0.8
lambda_02 <- 0.5
lambda_03 <- 0.3
lambda_04 <- 0.3
lambda_05 <- 0.2

pars <- c(n, mu, var_e, 
          beta_01, lambda_01, beta_02, lambda_02, beta_03, lambda_03,
          beta_04, lambda_04, beta_05, lambda_05)

Month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

          
### Generate simulation data
simulate_y <- function(pars, seed = 123) {
   set.seed(seed)
   
   ## Simulation parameters
   n         <- pars[1] # num of observation
   mu        <- pars[2] # intercept
   var_e     <- pars[3] # residual variance
   
   beta_01   <- pars[4] # regression coefficient of X1 to be esitmated
   lambda_01 <- pars[5] # decay rate of Ad-Stock effect of X1
   beta_02   <- pars[6]
   lambda_02 <- pars[7]
   beta_03   <- pars[8]
   lambda_03 <- pars[9]
   beta_04   <- pars[10]
   lambda_04 <- pars[11]
   beta_05   <- pars[12]
   lambda_05 <- pars[13]
   
   ## Create true Ad-Stock variables
   X_01_raw <- rgamma(n, 3) * ifelse(runif(n) > 0.7, 0, 1)
   X_01_fil <- stats::filter(X_01_raw, lambda_01, "recursive")
   
   X_02_raw <- rgamma(n, 2) * ifelse(runif(n) > 0.8, 0, 1)
   X_02_fil <- stats::filter(X_02_raw, lambda_02, "recursive")

   X_03_raw <- rgamma(n, 1) * ifelse(runif(n) > 0.4, 0, 1)
   X_03_fil <- stats::filter(X_03_raw, lambda_03, "recursive")

   X_04_raw <- rgamma(n, 1) * ifelse(runif(n) > 0.2, 0, 1)
   X_04_fil <- stats::filter(X_04_raw, lambda_04, "recursive")

   X_05_raw <- rgamma(n, 1) * ifelse(runif(n) > 0.5, 0, 1)
   X_05_fil <- stats::filter(X_05_raw, lambda_05, "recursive")
   
   ## Create Year & Seasonality
   year_eff_all     <- rep(runif(n_year) , each = n_per_year)
   seasonal_eff     <- rep(sin(-2:9), each = 4)
   seasonal_eff_all <- rep(seasonal_eff, n_year)

   ## Create residuals
   error <- rnorm(n, 0, sqrt(var_e))
   
   ## Create observations   
   y     <- mu + year_eff_all + seasonal_eff_all + 
      beta_01 * X_01_fil + beta_02 * X_02_fil + beta_03 * X_03_fil + 
      beta_04 * X_04_fil + beta_05 * X_05_fil + error
   
   ## Return dataset
   dat <- data.frame(
      "Y"          = y,
      "X_01"       = X_01_raw,
      "X_02"       = X_02_raw,
      "X_03"       = X_03_raw,
      "X_04"       = X_04_raw,
      "X_05"       = X_05_raw,
      "X_01_Fil"   = X_01_fil,
      "X_02_Fil"   = X_02_fil,
      "X_03_Fil"   = X_03_fil,
      "X_04_Fil"   = X_04_fil,
      "X_05_Fil"   = X_05_fil,
      "Y_lag"      = dplyr::lag(y, 1),
      "True_Error" = error,
      "Year"       = rep(2001:(2001 + n_year - 1), each = n_per_year * n_per_month),
      "Month"      = rep(rep(Month, each = n_per_month), n_year))
   return(dat)
}

## Data simulation
dat  <- na.omit(simulate_y(pars))


### reg for GA ##############################
my_regression <- function(dat, pred_vars) {
   tmp       <- dat[, c("Y", pred_vars)]
   model     <- glm(Y ~ ., tmp, family = "gaussian")
   return(model)
}

### mod and ROI function ##############################
get_model <- function(x) {
   index_selected     <- which(x == 1, arr.ind = T)
   pred_vars_selected <- pred_vars_all[index_selected]
   pred_vars          <- c(fixed_vars, pred_vars_selected)
   
   mod <- my_regression(dat, pred_vars)
   return(1.0 / mod$aic)
}

### AIC function ##############################
# get_AIC <- function(x) {
#    return(1.0 / mod$aic)
# }

### fixed variables ##############################
fixed_vars <- list(
   "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
   "FiscalYearEnd", "GW", "OBON", "SW", "nenmatsunenshi",
   "NEWYEAR.DAYS", "HolidayNum", "TaxRateUp"
)

### variables to be optimized ##############################
pred_vars_all <- list(
   "TV", "Digital", "Flyer", "Newspaper", "Magazine", "Radio",
   "Campaign_01", "Campaign_02", "Campaign_03"
)

### initial settings for genetic algorithm
population_size <- 20
generation_num  <- 20

### elite solution (best found so far)
elite_solution <- c(1, 1, 1, 0, 0, 0, 0, 1, 0)

### similar solutions to eliteSolution
other_solution <- matrix(elite_solution, population_size - 1, length(pred_vars_all),
                         byrow = T)
for (i in 1:(population_size - 1)) {
   for (j in 1:length(pred_vars_all)) {
      if (runif(1) < 1.0 / length(pred_vars_all)) {
         if (other_solution[i, j] == 0) other_solution[i, j] <- 1
         else otherSolution[i, j] <- 0
      }
   }
}

### initial population is generated with eliteSolution and otherSolution
# init_pop <- rbind(eliteSolution, otherSolution)

### genetic operations and parameters are set appropriately
gaControl("binary" = list(selection = "gabin_tourSelection",
                          crossover = "gabin_uCrossover"))

### genetic algorithm running
run_GA <- ga(type = "binary",
             fitness = get_model,
             # suggestions = initPop,
             pmutation = ga_pmutation,
             popSize = population_size,
             maxiter = generation_num,
             nBits = length(pred_vars_all),
             monitor = T)

### summary shows up and the data is saved
summary(runGA)
opt_indexselected <- which(runGA@solution[1,] == 1, arr.ind = T)
opt_predVarsSelected <- predVarsAll[opt_indexselected]
opt_fitness <- runGA@fitnessValue
opt_aic <- 1.0 / opt_fitness
opt_ModRoi <- getModRoi(runGA@solution[1,])
opt_result <- list(opt_predVarsSelected, opt_aic, opt_ModRoi[1][[1]], opt_ModRoi[2][[1]])
save(opt_result, file="GA_SizePopXXX_NumGenXXX.Rdata")
