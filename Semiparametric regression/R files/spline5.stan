data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, 2] X;
  matrix[N, K] Z;
}
parameters {
  vector[2] beta;
  vector[K] mspl;
  real<lower=0> sigma;
  real<lower=0> sigmaspl;
}
model {
  for (i in 1:N) {
   y[i] ~ normal(X[i,] * beta[1:2] + Z[i,] * mspl[1:K],sigma);   
}
for (k in 1:K){
  mspl[k] ~ normal(0,sigmaspl);
}
for (j in 1:2) {
  beta[j] ~ normal( 0  ,1000);
}
}

