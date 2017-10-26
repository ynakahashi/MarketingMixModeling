################################################################################
##
## Ad-stock estimation by rethinking
## written by Y.Nakahashi 
## 2017/10/26
##
################################################################################

################################################################################
## Environmental Settings
################################################################################
## Set working directory
# work_dir <- "/Users/nakahashi/Desktop/GitTest/MarketingMixModeling"
# setwd(work_dir)

## Load libraries
# devtools::install_github("rmcelreath/rethinking")
library(rethinking)

## Define functions
Create_Filtered_Vars <- function(pars, dat = Dat_Ori) {
   
   Dat_Tmp     <- dat
   len         <- length(pars)
   Target_Vars <- names(pars)
   
   for (i in 1:len) {
      Var <- Target_Vars[i]
      Fil <- stats::filter(Dat_Tmp[, Var], pars[i], "recursive")
      Dat_Tmp[Var] <- Fil
   }
   
   return(Dat_Tmp)
}

################################################################################
## Generate ample data
################################################################################
## Initial values
set.seed(456)
X1 <- runif(100); X2 <- runif(100); X3 <- runif(100); X4 <- runif(100); X5 <- runif(100);
l1 <- 0.8; l2 <- 0.6; l3 <- 0.5; l4 <- 0.8; l5 <- 0.7
b1 <- 0.5; b2 <- 0.8; b3 <- 0.3; b4 <- 0.2; b5 <- 0.6
bs <- c(3, b1, b2, b3, b4, b5)

Dat_Try             <- data.frame(cbind(X1, X2, X3, X4, X5))
Candidate_Vars      <- colnames(Dat_Try)
Lambda_Table        <- c(l1, l2, l3, l4, l5)
names(Lambda_Table) <- colnames(Dat_Try)

Dat_Fil   <- Create_Filtered_Vars(Lambda_Table, dat = Dat_Try)
Dat_Fil$Y <- Dat_Try$Y <- as.matrix(cbind(1, Dat_Fil)) %*% bs + rnorm(100, 0, 0.2)

## Check regression coefficients
summary(lm(Y ~ ., Dat_Fil))

################################################################################
## Estimate
################################################################################













