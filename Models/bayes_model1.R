rm(list=ls())

# Data Directory
#setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
# save.image(file = "bayes_m1yr.RData")
# load("bayes_m1yr.RData")

#Libraries
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)

# Load Data
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)
# Standardise - centre and scale data
X[,2:5] <- scale(X[,2:5], center = T, scale=T)
# Format datetime
X$datetime <- as.POSIXct(X$datetime, tz="GMT")
y_nov$datetime <- as.POSIXct(y_nov$datetime, tz="GMT")
y_pet$datetime <- as.POSIXct(y_pet$datetime, tz="GMT")
y_row$datetime <- as.POSIXct(y_row$datetime, tz="GMT")
# Merge predictors and variables into data frame
dat_nov <- merge(X, y_nov, by='datetime')
dat_nov <- dat_nov[,1:6]
dat_pet <- merge(X, y_pet, by='datetime')
dat_pet <- dat_pet[,1:6]
dat_row <- merge(X, y_row, by='datetime')
dat_row <- dat_row[,1:6]
colnames(dat_nov) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resNov')
colnames(dat_pet) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resPet')
colnames(dat_row) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resRow')

dim(dat_nov)


# Define Stan model
stan_code <- "
data {
  int<lower=0> N; // number of observations
  vector[N] response; // response variable
  matrix[N, 4] predictors; // time series predictors
}

parameters {
  real intercept;
  vector[4] beta;
  real<lower=0> sigma; // standard deviation of the error
}

model {
  intercept ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  
  for (i in 1:N) {
    response[i] ~ normal(intercept + dot_product(predictors[i], beta), sigma);
  }
}

generated quantities {
  vector[N] predicted;
  for (i in 1:N) {
    predicted[i] = normal_rng(intercept + dot_product(predictors[i], beta), sigma);
  }
}
"

# Compile the Stan model
stan_model <- cmdstan_model(write_stan_file(stan_code))

# Split Training Data set
train_dat <- dat_nov[1:8760,] ## 26204
train_dt <- train_dat[,1]
train_resp <- as.vector(train_dat[,6])
train_dat <- as.matrix(train_dat[,2:5])

# Prepare data for Stan model
stan_data <- list(
  N = nrow(train_dat),
  response = train_resp,
  predictors = train_dat # exclude response variable
)

# Fit the model to the data
fit <- stan_model$sample(data = stan_data)

# Print summary of the fit
print(fit)

# Plot posterior distributions
library(bayesplot)
posterior_samples <- fit$draws()
mcmc_pairs(posterior_samples, pars = c("intercept", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "sigma"))




