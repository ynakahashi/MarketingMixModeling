### library ##############################
library(GA)

### Parameters settings
n_year  <- 4
n_month <- 12
n_week  <- 4

n       <- n_year * n_month * n_week
mu      <- 4
var_e   <- 0.5

b_TV    <- 0.50
b_DG    <- 0.35 # Digital
b_Flyer <- 0.20
b_NP    <- 0.00 # Newspaper
b_MG    <- 0.00 # Magazine
b_Radio <- 0.00

l_TV    <- 0.8
l_DG    <- 0.5
l_Flyer <- 0.3
l_NP    <- 0.3
l_MG    <- 0.2
l_Radio <- 0.1

b_Camp_01 <- 0.3
b_Camp_02 <- 0.0
b_Camp_03 <- 0.4

pars <- c(n, mu, var_e, 
          b_TV, l_TV, b_DG, l_DG, b_Flyer, l_Flyer, b_NP, l_NP, b_MG, l_MG,
          b_Radio, l_Radio, b_Camp_01, b_Camp_02, b_Camp_03)

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
   beta_06   <- pars[14]
   lambda_06 <- pars[15]
   
   beta_11   <- pars[16]
   beta_12   <- pars[17]
   beta_13   <- pars[18]
   
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

   X_06_raw <- rgamma(n, 1) * ifelse(runif(n) > 0.2, 0, 1)
   X_06_fil <- stats::filter(X_06_raw, lambda_06, "recursive")
   
   ## Create Campaign effect
   X_11_raw <- ifelse(runif(n) > 0.2, 0, 1)
   X_12_raw <- ifelse(runif(n) > 0.6, 0, 1)
   X_13_raw <- ifelse(runif(n) > 0.4, 0, 1)
   
   ## Create Year & Seasonality
   year_eff_all     <- rep(runif(n_year) , each = n_month * n_week)
   seasonal_eff     <- rep(sin(-2:9), each = n_week)
   seasonal_eff_all <- rep(seasonal_eff, n_year)

   ## Create residuals
   error <- rnorm(n, 0, sqrt(var_e))
   
   ## Create observations   
   y     <- mu + year_eff_all + seasonal_eff_all + 
      beta_01 * X_01_raw + beta_02 * X_02_raw + beta_03 * X_03_raw + 
      beta_04 * X_04_raw + beta_05 * X_05_raw + beta_06 * X_06_raw + 
      beta_11 * X_11_raw + beta_12 * X_12_raw + beta_13 * X_13_raw + error
   
   # y     <- mu + year_eff_all + seasonal_eff_all + 
   #    beta_01 * X_01_fil + beta_02 * X_02_fil + beta_03 * X_03_fil + 
   #    beta_04 * X_04_fil + beta_05 * X_05_fil + beta_06 * X_06_fil + 
   #    beta_11 * X_11_raw + beta_12 * X_12_raw + beta_13 * X_13_raw + error
   
   ## Return dataset
   dat <- data.frame(
      "Y"              = y,
      "Year"           = paste0("Year_", rep(2001:(2001 + n_year - 1), each = n_month * n_week)),
      "Month"          = factor(rep(rep(Month, each = n_week), n_year),
                                levels = Month),
      "TV"             = X_01_raw,
      "Digital"        = X_02_raw,
      "Flyer"          = X_03_raw,
      "Newspaper"      = X_04_raw,
      "Magazine"       = X_05_raw,
      "Radio"          = X_06_raw,
      "Campaign_01"    = X_11_raw,
      "Campaign_02"    = X_12_raw,
      "Campaign_03"    = X_13_raw,
      "TV_Fil"         = X_01_fil,
      "Digital_Fil"    = X_02_fil,
      "Flyer_Fil"      = X_03_fil,
      "Newspaper_Fil"  = X_04_fil,
      "Magazine_Fil"   = X_05_fil,
      "Radio_Fil"      = X_06_fil,
      "Y_lag"          = dplyr::lag(y, 1),
      "True_Year_Eff"  = year_eff_all,
      "True_Month_Eff" = seasonal_eff_all,
      "True_Error"     = error)
   return(dat)
}

## Data simulation
dat  <- na.omit(simulate_y(pars))
dat_Year  <- psych::dummy.code(dat$Year)
dat_Month <- psych::dummy.code(dat$Month)

dat_dummy <- data.frame(cbind(
   dat[, -which(colnames(dat) %in% c("Year", "Month"))],
   dat_Year,
   dat_Month))


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
   
   mod <- my_regression(dat_Ana, pred_vars)
   return(1.0 / mod$aic)
}

### fixed variables ##############################
fixed_vars <- c(
   "Year_2001", "Year_2002", "Year_2003", "Year_2004", 
   "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

### variables to be optimized ##############################
pred_vars_all <- c(
   "TV", "Digital", "Flyer", "Newspaper", "Magazine", "Radio",
   "Campaign_01", "Campaign_02", "Campaign_03"
)

### data for analysis
dat_Ana <- dat_dummy[, c("Y", fixed_vars, pred_vars_all)]

### initial settings for genetic algorithm
population_size <- 50
generation_num  <- 200

### elite solution (best found so far)
elite_solution <- c(1, 1, 1, 0, 0, 0, 1, 0, 1)

### similar solutions to eliteSolution
other_solution <- matrix(elite_solution, population_size - 1, length(pred_vars_all),
                         byrow = T)


set.seed(123)
NS  <- nrow(other_solution) * ncol(other_solution)
r   <- population_size - 1
c   <- length(pred_vars_all)
idx <- matrix(runif(NS), r, c, byrow = T) < (1.0 / c)
x <- other_solution[idx]
other_solution[idx] <- as.integer(xor(x, TRUE))

### initial population is generated with eliteSolution and otherSolution
init_pop <- rbind(eliteSolution, otherSolution)

### genetic operations and parameters are set appropriately
gaControl("binary" = list(selection = "gabin_tourSelection",
                          crossover = "gabin_uCrossover"))

### genetic algorithm running
run_GA <- ga(type        = "binary",
             fitness     = get_model,
             suggestions = init_pop,
             pmutation   = ga_pmutation,
             popSize     = population_size,
             maxiter     = generation_num,
             nBits       = length(pred_vars_all),
             monitor     = T)

sum(run_GA@solution == elite_solution)

dat_Tmp <- dat_Ana[, c("Y", fixed_vars, pred_vars_all[run_GA@solution == 1])]
res_Tmp <- lm(Y ~ ., dat_Tmp)
coef(res_Tmp)
plot(coef(res_Tmp)[2:4], 
     (unique(dat$True_Year_Eff) - unique(dat$True_Year_Eff)[4])[-4])
plot(coef(res_Tmp)[6:16], 
     (unique(dat$True_Month_Eff) - unique(dat$True_Month_Eff)[12])[-12])



