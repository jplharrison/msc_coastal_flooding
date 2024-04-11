data {
  int<lower=0> N;
  array[N] real<lower=0,upper=1> NovLevel;
  array[N] real<lower=990,upper=1040> AtmPres;
  array[N] real<lower=-11,upper=14> LocalWind;
  array[N] real<lower=0,upper=17000> NapaFlow;
  array[N] real<lower=-14,upper=16> OceanWind;
}
parameters {
  real<lower=0,upper=1> beta0;
  real<lower=0,upper=1> beta1;
  real<lower=0,upper=1> beta2;
  real<lower=0,upper=1> beta3;
  real<lower=0,upper=1> beta4;
  real<lower=0,upper=1> sigma2;
}
model {
  theta ~ beta(1,1);  // uniform prior on interval 0,1
  y ~ bernoulli(theta);
}
