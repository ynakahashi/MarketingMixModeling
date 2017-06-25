data {
   int N;
   // int K = 1;
   vector[N] Y;
   vector[N] X_01;
   vector[N] X_02;
}

parameters {
   // real state_t0;
   real mu;
   real beta_X_01;
   real beta_X_02;
   real<lower = 0, upper = 1> lambda_X_01;
   real<lower = 0, upper = 1> lambda_X_02;
   // vector[N*K] state_t;
   // vector[N*K] beta_X;
   // real<lower=0> var_state_t0;
   // real<lower=0> var_state;
   // real<lower=0> var_beta_X;
   real<lower=0> var_error;
   real<lower=0> var_lambda_01;
   real<lower=0> var_labmda_02;
}

transformed parameters {
   real X_01_fil[N];
   real X_02_fil[N];
   
   ## Create Ad-stock variable
   for(j in 1:(N)) {
      if(j == 1) {
         X_01_fil[j] = X_01[j];
         X_02_fil[j] = X_02[j];
      } else {
         // X_01_fil[j] = X_01[j] + X_01_fil[j-1] ^ lambda_X_01;
         // X_02_fil[j] = X_02[j] + X_02_fil[j-1] ^ lambda_X_02;
         X_01_fil[j] ~ normal(X_01[j] + X_01_fil[j-1] ^ lambda_X_01[j],
         var_lambda_01);
         X_02_fil[j] ~ normal(X_02[j] + X_02_fil[j-1] ^ lambda_X_02[j],
         var_labmda_02);
      }
   }
}

model {

   ## Order Num
   for(i in 1:(N)) {
      Y[i] ~ normal(mu + beta_X_01 * X_01_fil[i] + beta_X_02 * X_02_fil[i], 
      var_error);
   }
}
