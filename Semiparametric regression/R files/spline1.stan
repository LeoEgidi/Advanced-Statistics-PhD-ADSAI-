data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> K;
  vector[N] y;
  matrix[N, 4] X;
  matrix[N, K] Z;
}
parameters {
  vector[4] beta;
  vector[K] b;
  real<lower=0> sigma2;
  real<lower=0> sigmab2;
}
model {
  for (i in 1:N) {
   y[i] ~ normal( X[i,]*beta[1:(p+1)] + Z[i,]*b[1:K],sqrt(sigma2));   
  }
  for (k in 1:K){
    b[k] ~ normal(0,sqrt(sigmab2));
  }
  for (j in 1:(p+1)) {
    beta[j] ~ normal( 0  ,1000);
  }
  sigma2 ~ inv_gamma(0.001,0.001);
  sigmab2 ~ inv_gamma(0.001, 0.001);
}
