rm(list=ls())

# Data Directory
#setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
# save.image(file = "bayes_stan.RData")
# load("bayes_stan.RData") #load("bayes_m1yr.RData")

#Libraries
library(bayesplot)
library(coda)
library(dplyr)
library(ggplot2)
library(rstan)
library(posterior)

# Load Data
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)

# Rescale predictors between 0 and 1
for(i in 2:5){X[,i] <- rescale(X[,i], to=c(0,1))}
# Format datetime
X$datetime <- as.POSIXct(X$datetime, tz="GMT")
y_nov$datetime <- as.POSIXct(y_nov$datetime, tz="GMT")
y_pet$datetime <- as.POSIXct(y_pet$datetime, tz="GMT")
y_row$datetime <- as.POSIXct(y_row$datetime, tz="GMT")

# Ocean_onshorewind: Ocean Wind, located at the NDBC buoy #46026,  On-shore component (from 100 N) 
# Gnoss_onshorewind: Local Wind at the Gnoss Field Airport, On-shore component (from 60 N) 
# AtmPres: Atmospheric Pressure, mBAR
# napa_flow_cfs: River flow at Napa Ricer at USGS gauge 11458000
# Field residual =  stage - predicted It represents the variation of the water level due to non-tidal forcing.
# Field stage_m: raw water level data - detrended by removing the mean value, and resampled to hourly time intervals. Small data gaps were filled by linear interpolation. 
# Field predicted_m:  This is the predicted tide. The predicted tide was calculated using a publicly available Python routine based on a well documented Matlab routine called Utide (http://www.po.gso.uri.edu/~codiga/utide/utide.htm). 

# Merge predictors and variables into data frame - RESIDUAL WATER LEVEL RESPONSE ~ 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resNov'
dat_nov <- merge(X, y_nov, by='datetime')
dat_nov <- dat_nov[,1:6]
dat_pet <- merge(X, y_pet, by='datetime')
dat_pet <- dat_pet[,1:6]
dat_row <- merge(X, y_row, by='datetime')
colnames(dat_nov) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resNov')
colnames(dat_pet) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resPet')
colnames(dat_row) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resRow')

# Create time series object
dt <- dat_nov[,1]
dat <- ts(data=dat_nov[2:6], start=c(2019, 20),
          end=c(2022,6477.562), frequency=365.24219*24, class="matrix")
dat <- as.data.frame(dat); 
# 32756     9
N <- dim(dat)[1]
N_test <- 2448
N_train <- N-N_test
dat_train <- dat[1:N_train,1:4]
resNov_train <- dat[1:N_train,5]
dt_train <- dt[1:N_train]
dat_test <- dat[(N_train+1):N,1:4]
resNov_test <- dat[(N_train+1):N,5]
dt_test <- dt[(N_train+1):N]
var_names <- c('Atmospheric Pressure', 'Local Wind', 'Napa Flow', 'Ocean Wind', 'Residual Water Level')


# Define Stan model
stan_code <- "
data {
  int<lower=0> N;                 // Number of observations
  int<lower=0> K;                 // Number of predictors
  matrix[N, K] X;                 // Predictor matrix
  vector[N] y;                    // Response variable
}

parameters {
  vector[K] beta;                 // Coefficients for predictors
  real alpha;                     // Intercept
  real<lower=0> sigma;            // Error standard deviation
}

model {
  y ~ normal(X * beta + alpha, sigma); // Likelihood
  alpha ~ normal(0, 10);               // Prior for intercept
  beta ~ normal(0, 1);                 // Prior for coefficients
  sigma ~ cauchy(0, 2);                // Prior for sigma
}
"

# data {
#   int<lower=0> N; // number of observations
#   vector[N] response; // response variable
#   matrix[N, 4] predictors; // time series predictors
# }
# 
# parameters {
#   real intercept;
#   vector[4] beta;
#   real<lower=0> sigma; // standard deviation of the error
# }
# 
# model {
#   intercept ~ normal(0, 10);
#   beta ~ normal(0, 10);
#   sigma ~ normal(0, 10);
#   
#   for (i in 1:N) {
#     response[i] ~ normal(intercept + dot_product(predictors[i], beta), sigma);
#   }
# }
# 
# generated quantities {
#   vector[N] predicted;
#   for (i in 1:N) {
#     predicted[i] = normal_rng(intercept + dot_product(predictors[i], beta), sigma);
#   }
# }
#stan_model <- cmdstan_model(write_stan_file(stan_code))
stan_model <- stan_model(model_code = stan_code)
## add seasonality, penalty (priors)

# Data for Stan model
stan_data <- list(
  N = N_train,
  K = dim(dat_train)[2],
  X = dat_train,
  y = resNov_train # exclude response variable
)

# Fit the model to the data
#fit <- stan_model$sample(data = stan_data)
fit <- sampling(
  stan_model,
  data = stan_data,
  iter = 2000,        # Number of iterations
  chains = 4,         # Number of chains
  warmup = 500,       # Number of warmup iterations
  thin = 1,           # Thinning interval
  seed = 123          # Random seed for reproducibility
)

#https://mc-stan.org/misc/warnings.html#bfmi-low 
#https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

# Print summary of the fit
print(fit)
summary(fit)

# Plot posterior distributions
library(bayesplot)
posterior_samples <- extract(fit)
beta_samples <- posterior_samples$beta
alpha_samples <- posterior_samples$alpha
sigma_samples <- posterior_samples$sigma


# Convert samples to a data frame
beta_df <- as.data.frame(posterior_samples$beta)
alpha_df <- as.data.frame(posterior_samples$alpha)
sigma_df <- as.data.frame(posterior_samples$sigma)

beta_samples <- as.matrix(fit, pars = c("alpha", "beta", "sigma"))
mcmc_list <- as.mcmc.list(lapply(1:ncol(beta_samples), function(i) {
  as.mcmc(beta_samples[, i])
}))
mcmc_pairs(mcmc_list, pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "alpha", "sigma"))
## want to see independence
mcmc_trace(as.array(fit), pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "alpha", "sigma"))

par(mfrow=c(3,2))
hist(posterior_samples$beta[,1])
hist(posterior_samples$beta[,2])
hist(posterior_samples$beta[,3])
hist(posterior_samples$beta[,4])

hist(posterior_samples$sigma)
hist(posterior_samples$alpha)

