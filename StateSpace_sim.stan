data {
   int N; // 地域ごとの観測値の数（data_length）
   vector[N] Y; // 観測値のベクトル
   vector[N] X; // 説明変数のベクトル
}

parameters {
   real state_t0;     // 状態変数の0期目の平均
   vector[N] state; // 状態変数のベクトル
   //   real beta_0;       // 回帰係数の事前分布の平均
   real beta;    // 地域ごとの回帰係数のベクトル
   real<lower=0> var_state_t0; // 状態変数の0期目の分散
   real<lower=0> var_state;    // 状態変数の分散
   //   real<lower=0> var_beta_0;   // 回帰係数の事前分布の分散
   real<lower=0> var_error;    // 誤差分散
}

model {
   // 状態変数をサンプリング
   // 1期目の値は0期目の分布からサンプルする
   state[1] ~ normal(state_t0, var_state_t0);
   
   // 2期目以降は前期の値を平均とした分布からサンプルする
   for(i in 2:N) {
      state[i] ~ normal(state[i-1], var_state);
   }

   // 回帰係数をサンプリング
   beta ~ normal(0, 1);
   
   // Yをサンプリング
   for(i in 1:N) {
      Y[i] ~ normal(state[i] + beta * X[i], var_error);
   }
}
