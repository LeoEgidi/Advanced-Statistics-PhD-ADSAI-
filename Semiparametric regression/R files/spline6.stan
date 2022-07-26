data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> KB;
  vector[N] y;
  int ind[N];
  matrix[N, 2] X;
  matrix[N, K] Z;
  matrix[N, KB] ZB;
}
parameters {
  vector[2] beta;
  matrix[2,2] suspbeta;
  matrix[2,KB] susp;
  vector[K] mspl;
  real<lower=0> sigma;
  real<lower=0> sigmaspl;
  real<lower=0> sigmau;
  real<lower=0> sigmasusp;
}
model {
  for (i in 1:N) {
   y[i] ~ normal(X[i,] * beta[1:2] + Z[i,] * mspl[1:K] + 
            X[i,] * suspbeta[ind[i],1:2]' + 
            ZB[i,] * susp[ind[i],1:KB]',sigma);   
}
for (k in 1:K){
  mspl[k] ~ normal(0,sigmaspl);
}
for (h in 1:2){ 
   suspbeta[h,1] ~ normal( 0  ,sigmau);
   suspbeta[h,2] ~ normal( 0  ,sigmau);  
   for (k in 1:KB){
    susp[h,k] ~ normal(0,sigmasusp);
   }}
for (j in 1:2) {
  beta[j] ~ normal( 0  ,1000);
}
}
