data {
   int N; // number of observation in each region
   int K; // number of region
   vector[N*K] Y; // vector for Y
   vector[N*K] X; // vector for X
   int<lower=1, upper=K> Area_ID[N*K]; // vector for Area ID
}

parameters {
   real state_t0;
   real beta_0;
   vector[N*K] state;
   vector[K] beta;
   real<lower=0> var_state_t0;
   real<lower=0> var_state;
   real<lower=0> var_beta_0;
   real<lower=0> var_error;
}


model {
   // state
   // for (k in 1:K) {
   //    state[1 + (k-1)*N] ~ normal(state_t0, var_state_t0);
   //    for(i in 2:N) {
   //       state[i + (k-1)*N] ~ normal(state[i-1 + (k-1)*N], var_state);
   //    }
   // }

   // stateの1期目は指定しない（岩波DS Vol.1の松浦さんのスクリプトを参考）
   // for (k in 1:K) {
   //    for(i in 2:N) {
   //       state[i + (k-1)*N] ~ normal(state[i-1 + (k-1)*N], var_state);
   //    }
   // }

   for (i in 1:N*K) {
      if (i % N == 1 ) {
         state[i] ~ normal(state_t0, var_state_t0);
      } else {
         state[i] ~ normal(state[i-1], var_state);
      }
   }

   // beta
   for (k in 1:K) {
      beta[k] ~ normal(beta_0, var_beta_0);
   }

   // Y
   for(i in 1:(N*K)) {
      Y[i] ~ normal(state[i] + beta[Area_ID[i]] * X[i], var_error);
   }
}
