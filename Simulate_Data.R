N <- 300
must_var <- 7
select_var <- 5
select_pos <- c(T, T, F, F, T)
lambda_table <- c(0.8, 0.3, 0.6, 0.5, 0.4)
lambda_pos <- c(T, T, F, F, F)


Simulate_Data(N, must_var, select_var, select_pos, lambda_table, lambda_pos,
              trend = T, seasonal = T, p_season = 50, seed = 1234)

Simulate_Data <- function(N,            # Number of obsevations
                          must_var,     # Number of must-included variables
                          select_var,   # Number of potential variables 
                          select_pos,   # Index determining whether included or not for potential vars(LOGICAL)
                          lambda_table, # Lagged-effect values
                          lambda_pos,   # Index determining whether having lagged-effect or not for potential vars(LOGICAL)
                          min_b  = 0.5, # Minimum of regression coefficient
                          max_b  = 1.5, # Maximum of regression coefficient
                          type = "additive", # Additive model       
                          error_var = 4, # Residual variance
                          trend = F,  # LOGICAL if include trend term
                          trend_type = "linear", # Type of trend term
                          trend_init = 5, # Initial value of trend term
                          trend_var  = 2, # Variance of white noise for trend term
                          seasonal = F,   # LOGICAL if include seasonal effect
                          p_season = 2,   # Period of a season
                          seed  = 123,    # Random seed
                          ...) {

   ## Error check
   # if (!is.integer(must_var) | !is.integer(select_var)) {
   #    stop("Input integer number for must var and select var")
   # }
   if (select_var != length(select_pos)) {
      stop("Number mismatch for select var & select pos")
   }
   if (length(lambda_pos) != length(lambda_table)) {
      stop("Number mismatch for lambda var & select pos")
   }
   if (select_var < length(lambda_table)) {
      stop("Number mismatch for select var & lambda var")
   }
   if (trend_type != "linear") {
      stop("Sorry, not implemented yet!")
   }
   if (p_season <= 1 | p_season >= N) {
      stop("Strange seasonal cycle input")
   }
   
   ## Simulate data
   var_num <- must_var + select_var
   d       <- data.frame(matrix(NA, N, var_num))

   ## Simulate trend term
   set.seed(seed)
   if (trend) {
      d$Trend <- NA
      d$Trend[1] <- trend_init + rnorm(1, 0, sqrt(trend_var))
      for (t in 2:N) {
         d$Trend[t] <- d$Trend[t - 1] + rnorm(1, 0, sqrt(trend_var))
      }
      var_trend <- "trend"
   } else {
      d$Trend <- 0
   }
   
   ## Simulate seasonal term
   if (seasonal) {
      d$season <- rep(runif(p_season), length.out = N)
      var_season <- "season"
   } else {
      d$season <- 0
   }
   
   ## Simulate X for fixed vars
   d[, 1:must_var] <- matrix(rnorm(N * must_var), N)
   var_must <- paste0("Must_", 1:must_var)
   
   ## Simulate X for select vars
   s_f <- s <- matrix(rgamma(N * select_var, 1), N)
   
   for (k in 1:select_var) {
      dr       <- lambda_pos[k] * lambda_table[k]
      s_f[, k] <- stats::filter(s_f[, k], dr, "recursive")
   }
   d[, (must_var + 1):var_num] <- s_f
   var_select <- paste0("select_", 1:select_var)
   
   ## Generate beta
   a     <- 0
   b_mst <- runif(must_var,   min = min_b, max = max_b)
   b_sel <- runif(select_var, min = min_b, max = max_b) * select_pos

   ## Generate Y
   d$Y <- a + d$Trend + d$season + as.matrix(d[, 1:var_num]) %*% c(b_mst, b_sel) + rnorm(N, 0, sqrt(error_var))
   colnames(d) <- c(var_must, var_select, "Trend", "season", "Y")
   
   ## Return Data
   Out <- list("Out_data" = d,
               "Parameters" = list(
                  "Intercept"       = a,
                  "Regression_Coef" =  c(b_mst, b_sel),
                  "Select_Position" = select_pos,
                  "Lambda_Value"    = lambda_table,
                  "Lambda_Position" = lambda_pos,
                  "Trend_Init"      = trend_init,
                  "Trend_Var"       = trend_var,
                  "Residual_Var"    = error_var,
                  "Seasonal_Cycle"  = p_season))
   return(Out)
}

