### library ##############################
library(GA)

### Parameters settings
n_year     <- 4
n_per_year <- 48 # 4 weeks * 12 months
n          <- n_year * n_per_year
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
   seasonal_eff     <- sin(seq(from = -3, to = 3, length.out = n_per_year))
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
      "True_Error" = error)
   return(dat)
}

## Data simulation
dat  <- na.omit(simulate_y(pars))


### reg for GA ##############################
regGA <- function(dat, pred_vars) {
   var_names <- pred_vars
   tmp       <- dat[, c("Y", var_names)]
   model     <- lm(Y ~ ., tmp)
   return(model)
}

### mod and ROI function ##############################
getModRoi <- function(x) {
   indexSelected    <- which(x == 1, arr.ind = T)
   predVarsSelected <- predVarsAll[indexSelected]
   predVars <- c(fixedVars, predVarsSelected)
   mod <- regGA(datLrn, predVars)
}


### AIC function ##############################
getFitness <- function(x) {
   return(1.0/mod$aic)
}

### fixed variables ##############################
fixedVars <- list(
   "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
   "Week_idx_02", "Week_idx_03", "Week_idx_04", "Week_idx_05", 
   "FiscalYearEnd", "GW", "OBON", "SW", "nenmatsunenshi",
   "NEWYEAR.DAYS", "HolidayNum", "TaxRateUp"
)

### variables to be optimized ##############################
predVarsAll <- list(
   "TV", "Digital", "Flyer", "Newspaper", "Magazine", "Radio",
   "Campaign_01", "Campaign_02", "Campaign_03"
)

### initial settings for genetic algorithm
Population_Size <- 20
Generation_Num  <- 20

### elite solution (best found so far)
eliteSolution <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                   0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 1)

### similar solutions to eliteSolution
otherSolution <- matrix(eliteSolution, populationSize - 1, length(predVarsAll), 
                        byrow = T)
for (i in 1:(populationSize - 1)) {
   for (j in 1:length(predVarsAll)) {
      if (runif(1) < 1.0 / length(predVarsAll)) {
         if (otherSolution[i, j] == 0) otherSolution[i, j] <- 1
         else otherSolution[i, j] <- 0
      }
   }
}

### initial population is generated with eliteSolution and otherSolution
initPop <- rbind(eliteSolution, otherSolution)

### genetic operations and parameters are set appropriately
gaControl("binary" = list(selection = "gabin_tourSelection",
                          crossover = "gabin_uCrossover"))

### genetic algorithm running
runGA <- ga(type = "binary",
            fitness = getFitness,
            suggestions = initPop,
            pmutation = ga_pmutation,
            popSize = populationSize,
            maxiter = generationNum,
            nBits = length(predVarsAll),
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