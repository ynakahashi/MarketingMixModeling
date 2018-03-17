data {
   int N;
   vector[N] Y;
   vector[N] X_01;
   vector[N] X_02;
}

parameters {
   real mu;
   real beta_01;
   real beta_02;
   real<lower=0, upper=1> lambda_01;
   real<lower=0, upper=1> lambda_02;
   real<lower=0, upper=1> k_01;
   real<lower=0, upper=1> k_02;
   real<lower=0, upper=5> s_01;
   real<lower=0, upper=5> s_02;
   real<lower=0> var_error;
}

transformed parameters {
   vector[N] X_01_fil;
   vector[N] X_02_fil;
   vector[N] X_01_conv;
   vector[N] X_02_conv;

   X_01_fil[1] = X_01[1];
   X_02_fil[1] = X_02[1];
   
   // Create Ad-stock variable
   for(j in 2:N) {
      X_01_fil[j] = X_01[j] + X_01_fil[j-1] * lambda_01;
      X_02_fil[j] = X_02[j] + X_02_fil[j-1] * lambda_02;
   }
   
   for(k in 1:N) {
      X_01_conv[k] = 1 / (1 + (X_01_fil[k] / k_01)^-s_01);
      X_02_conv[k] = 1 / (1 + (X_02_fil[k] / k_02)^-s_02);
   }
}

model {
   // Sampling
   for(i in 1:N) {
      Y[i] ~ normal(mu + beta_01 * X_01_conv[i] + beta_02 * X_02_conv[i], var_error);
   }
}
