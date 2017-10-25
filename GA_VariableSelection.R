################################################################################
##
## Variable selection & Ad-stock Estimation
## written by Y.Nakahashi 
## 2017/10/25
##
################################################################################

################################################################################
## Environmental Settings
################################################################################
## Set working directory
# work_dir <- "/Users/nakahashi/Desktop/GitTest/MarketingMixModeling"
# setwd(work_dir)

## Load libraries
library(GA)


################################################################################
## Create functions
################################################################################
## Update lambda values(Ad-stock rate)
## Input  : Variable Index (Indicates Include / Exclude)
## Output : Updated Lambda table under given Variable Index
Update_Lambda_Table <- function(Var_Idx, opt_method) {

   ## Select target variable & their lambda  
   Vars    <- Candidate_Vars[Var_Idx == 1]
   Lambda  <- Lambda_Table[which(names(Lambda_Table) %in% Vars)]
   
   ## Optimize lamabda
   if (opt_method == "optim") {
      
      ## by optim
      Values  <- optim(par = Lambda, 
                       fn = Return_AIC, 
                       method = "L-BFGS-B", 
                       lower = rep(0, length(Lambda)),
                       upper = rep(1, length(Lambda)))$par
      Values  <- optim(par = Values, 
                       fn = Return_AIC,
                       method = "L-BFGS-B", 
                       lower = rep(0, length(Lambda)),
                       upper = rep(1, length(Lambda)))
      
      ## Return updated lambda
      Res     <- Values$par 
      
   } else {
      
      ## by GA
      Population_Size <- 20
      Generation_Num  <- 100
      Values  <- ga(type = 'real-valued', 
                    min = rep(0, length(Lambda)), 
                    max = rep(1, length(Lambda)),
                    popSize = Population_Size, 
                    maxiter = Generation_Num, 
                    pmutation   = ga_pmutation,
                    names = Vars,
                    keepBest = T, 
                    monitor  = F,
                    fitness = Return_AIC)
      
      ## Return updated lambda
      Res     <- apply(Values@solution, 2, median)
      
   }
   
   ## Return updated lambda
   Lambda_Up <- Lambda_Table
   Lambda_Up[which(names(Lambda_Up) %in% Vars)] <- Res
   return(Lambda_Up)
}


## Return AIC under given lambda (Objective function)
## Input  : Lambda values
## Output : AIC value
Return_AIC <- function(pars) {
   Dat_Tmp <- Create_Filtered_Vars(pars, Dat_Ori)
   mod     <- glm(Y ~ ., data = Dat_Tmp, family = "gaussian")
   return(-mod$aic)
}

## Create filter-ed variable
## Input  : Lambda values
## Output : Data with filter-ed variables by given lambda
Create_Filtered_Vars <- function(pars, dat = Dat_Ori) {

   Dat_Tmp     <- dat
   len         <- length(pars)
   Target_Vars <- names(pars)

   for (i in 1:len) {
      Var <- Target_Vars[i]
      # if (is.na(Var)) next
      Fil <- stats::filter(Dat_Tmp[, Var], pars[i], "recursive")
      Dat_Tmp[Var] <- Fil
   }
   
   return(Dat_Tmp)
}


## Return Variable Index
## Input  : Lambda Table
## Output : Variable Index
Update_Variable_Index <- function(Lambda_Table, Pop_Size = 20, Gen_N = 100) {

   ### Create filter-ed variable
   Dat_Tmp <- Create_Filtered_Vars(Lambda_Table)
   
   ### Set options for GA
   Population_Size <- Pop_Size
   Generation_Num  <- Gen_N
   
   ### elite solution (best found so far)
   Elite_Sol <- Var_Idx
   
   ### similar solutions to eliteSolution
   Other_Sol <- matrix(Elite_Sol, Population_Size - 1, length(Candidate_Vars),
                       byrow = T)
   
   # set.seed(123)
   NS  <- prod(dim(Other_Sol))
   r   <- Population_Size - 1
   c   <- length(Candidate_Vars)
   idx <- matrix(runif(NS), r, c, byrow = T) < (1.0 / c)
   x   <- Other_Sol[idx]
   Other_Sol[idx] <- as.integer(xor(x, TRUE))

   if (any(rowSums(Other_Sol) == 0)) {
      Other_Sol[rowSums(Other_Sol) == 0, ] <- as.integer(runif(c) < 0.8)
   }
      
   ### initial population is generated with eliteSolution and otherSolution
   Init_Pop <- rbind(Elite_Sol, Other_Sol)
   
   ### genetic operations and parameters are set appropriately
   gaControl("binary" = list(selection = "gabin_tourSelection",
                             crossover = "gabin_uCrossover"))
   
   ### genetic algorithm running
   Variable_Index <- ga(type        = "binary",
                        fitness     = Return_AIC_GA,
                        suggestions = Init_Pop,
                        pmutation   = ga_pmutation,
                        popSize     = Population_Size,
                        maxiter     = Generation_Num,
                        nBits       = length(Candidate_Vars),
                        monitor     = T,
                        keepBest    = T,
                        Lambda_Table = Lambda_Table)
   
   return(Variable_Index@solution)
}


## Return AIC under given Variable Index & lambda (Objective function)
## Input  : Variable Index, Lambda values
## Output : AIC value
Return_AIC_GA <- function(x, Lambda_Table = Lambda_Table) {

   Vars      <- Candidate_Vars[x == 1]
   # Pred_Vars          <- c(Fixed_Vars, Vars)
   Pred_Vars <- c(Vars)
   Lambda    <- Lambda_Table[which(names(Lambda_Table) %in% Pred_Vars)]
   
   Return_AIC(Lambda)
}


myprintf <- function(fmt, ...) {
   cat(sprintf(fmt, ...))
}

################################################################################
## Try using sample data
################################################################################
## Initial values
set.seed(456)
X1 <- runif(100); X2 <- runif(100); X3 <- runif(100); X4 <- runif(100); X5 <- runif(100);
l1 <- 0.4; l2 <- 0.3; l3 <- 0.2; l4 <- 0.5; l5 <- 0.7
b1 <- 0.5; b2 <- 0.8; b3 <- 0.3; b4 <- 0.2; b5 <- 0.6
bs <- c(3, b1, b2, b3, b4, b5)

Dat_Try             <- data.frame(cbind(X1, X2, X3, X4, X5))
Candidate_Vars      <- colnames(Dat_Try)
Lambda_Table        <- c(l1, l2, l3, l4, l5)
names(Lambda_Table) <- colnames(Dat_Try)

Dat_Fil   <- Create_Filtered_Vars(Lambda_Table, dat = Dat_Try)
Dat_Try$Y <- as.matrix(cbind(1, Dat_Fil)) %*% bs + rnorm(100)
Dat_Ori   <- Dat_Try


## 
Lambda_Table        <- rep(0.2, length(Candidate_Vars))
names(Lambda_Table) <- colnames(Dat_Try)[1:length(Candidate_Vars)]
Var_Idx             <- rep(1, length(Candidate_Vars))
AIC_Old             <- -200
for (i in 1:3) {
   # cat("\n")
   AIC_Old      <- AIC_New
   myprintf("Update Lambda Table: %i", i)
   Lambda_Table <- Update_Lambda_Table(Var_Idx, opt_method = "optim")
   myprintf(" ... Complete \n")
   myprintf("Update Variable Index: %i", i)
   Res_GA       <- Update_Variable_Index(Lambda_Table)
   myprintf(" ... Complete \n")
   Var_Idx      <- Res_GA@solution
   AIC_New      <- Res_GA@fitnessValue
   myprintf("Lambda Table: %f \n", Lambda_Table)
   myprintf("Selected Variable: %i \n", Var_Idx)
   myprintf("\n AIC_New - AIC_Old: %f \n", (AIC_New - AIC_Old))
}




################################################################################
## Try using simulated data
################################################################################



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
plot((unique(dat$True_Year_Eff) - unique(dat$True_Year_Eff)[4])[-4],
     coef(res_Tmp)[2:4])
abline(a = 0, b = 1, col = "red", lty = 2)
plot((unique(dat$True_Month_Eff) - unique(dat$True_Month_Eff)[12])[-12],
     coef(res_Tmp)[6:16])
abline(a = 0, b = 1, col = "red", lty = 2)
plot(c(b_TV, b_DG, b_Flyer, b_Camp_01, b_Camp_02, b_Camp_03),
     coef(res_Tmp)[18:23])
abline(a = 0, b = 1, col = "red", lty = 2)
