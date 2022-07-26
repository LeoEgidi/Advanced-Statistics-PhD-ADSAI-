data {
  int<lower=0> N;
  int<lower=0> K;
  int y[N];
  matrix[N, 2] X;
  matrix[N, K] Z;
}
parameters {
  vector[2] beta;
  vector[K] b;
  real<lower=0> sigmab;
}
model {
  for (i in 1:N) {
   y[i] ~ bernoulli_logit(X[i,] * beta[1:2] + Z[i,] * b[1:K]);
}
for (k in 1:K){
  b[k] ~ normal(0,sigmab);
}
sigmab ~ normal(0,1000); 
for (j in 1:2) {
  beta[j] ~ normal( 0  ,1000);
}
}

