data {
  int<lower=0> N;
  array[N] real<lower=0,upper=1> NovLevel;
  array[N] real<lower=0,upper=1> AtmPres;
  array[N] int<lower=0,upper=1> LocalWind;
  array[N] int<lower=0,upper=1> NapaFlow;
  array[N] int<lower=0,upper=1> OceanWind;
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  theta ~ beta(1,1);  // uniform prior on interval 0,1
  y ~ bernoulli(theta);
}
