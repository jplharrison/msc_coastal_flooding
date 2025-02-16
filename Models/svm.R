rm(list=ls())


# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='svm.RData')
#save.image(file = 'svm-model.RData')
#save.image(file='svm-nu-tuning.RData')
#load(file = 'svm-model.RData')
#load(file='svm-nu-tuning.RData')

#save.image(file='svm-3.RData')
#load(file='svm-3.RData')

# save.image(file='svm-24-0-3.RData')

#Libraries
library(e1071)
library(forecast)
library(ggplot2)
library(gridExtra)
library(mlr)
library(scales)
library(tidyverse)
library(tseries)

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
dat_nov <- dat_nov[,1:8]
colnames(dat_nov) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resNov', 'stageNov', 'predNov')
dat_nov <- dat_nov %>% dplyr::select(datetime, atmPres, gnossWind, napaFlow, oceanWind, predNov, stageNov)
# Generate lagged features
N <- dim(X)[1]; N_test <- 2448; N_train <- N-N_test; N_lags=24
dat_lags <- matrix(NA, nrow=N, ncol=N_lags)
for(i in 1:N_lags){dat_lags[(1+i):(N),i]=dat_nov$stageNov[(1):(N-i)]}
colnames(dat_lags) <- paste0('lag_',1:N_lags)

# Test and Training Split
dat_train <- na.omit(cbind(dat_lags[1:N_train,], dat_nov[1:N_train,]))
dat_test <- cbind(dat_lags[(1+N_train):N,], dat_nov[(1+N_train):N,])
dt_train <- dat_nov$datetime[(1+N_lags):N_train]
dt_test <- dat_nov$datetime[(1+N_train):N]
stage_train <- dat_train$stageNov
stage_test <- dat_test$stageNov
dat_train <- as.data.frame(ts(data=(dat_train%>%select(colnames(dat_lags), atmPres, gnossWind, napaFlow, oceanWind, predNov)), start=c(2019, 20+N_lags),
                        end=c(2022,(6477.562-N_test)), frequency=365.24219*24, class="matrix"))
dat_test <- as.data.frame(na.omit(ts(data=(dat_test%>%select(colnames(dat_lags), atmPres, gnossWind, napaFlow, oceanWind, predNov)), start=c(2022,(6477.562-N_test+1)),
                       end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix")))
var_names <- c('Atmospheric Pressure', 'Local Wind', 'Napa Flow', 'Ocean Wind', 'Predicted Tide', 'Water Level')

## SVR24 Parameter Tuning
set.seed(890)
dat_train$stage <- stage_train
task <- makeRegrTask(data=dat_train, target='stage')
learner <- makeLearner('regr.svm', type='nu-regression')
param_set <- makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),#, "linear", "polynomial")),
  makeDiscreteParam("cost", values = c(2.5)),
  makeDiscreteParam("nu", values = c(0.5, 0.6)),
  makeDiscreteParam("gamma", values = c(0, 0.1))
)
tune_ctrl <- makeTuneControlGrid()
print(Sys.time())
tuned_result <- tuneParams(learner, task, resampling = makeResampleDesc("CV", iters = 5),
                           par.set = param_set, control = tune_ctrl)
print(Sys.time())

# SVR24 Model Fit
set.seed(4321)
lag_terms <- paste(paste0('lag_',1:24), collapse = " + ")
svm_formula <- as.formula(paste("stage_train ~ atmPres + gnossWind + napaFlow + oceanWind + predNov +", lag_terms))

svm_models <- list(); svm_rmses <- vector(); svm_rmsesIS <- vector()
i=1; print(Sys.time())
svm_models[[i]] <- svm(formula=svm_formula, 
                 data = dat_train, 
                 type= 'nu-regression',
                 kernel = "radial", # 'linear', 'sigmoid', 'polynomial', 'rbf' (radial)
                 nu = 0.6,
                 cost = 2.5, 
                 gamma = 0.1,
                 cross=0)
print(Sys.time())


set.seed(890)
print(Sys.time())
dat_train_3 <- dat_train %>% select(lag_1, lag_2, lag_3, atmPres, gnossWind, napaFlow , oceanWind, predNov)
dat_train_3$stage <- stage_train
task <- makeRegrTask(data=dat_train_3, target='stage')
learner <- makeLearner('regr.svm', type='nu-regression')
param_set <- makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),#, "linear", "polynomial")),
  makeDiscreteParam("cost", values = c(10)),
  makeDiscreteParam("nu", values = c(0.5, 0.6)),
  makeDiscreteParam("gamma", values = c(0, 0.1))
  # makeDiscreteParam("kernel", values = c("radial")),#, "linear", "polynomial")),
  # makeDiscreteParam("cost", values = c(1, 2.5, 5)),
  # makeDiscreteParam("nu", values = c(0.1, 0.3, 0.5)),
  # makeDiscreteParam("gamma", values = c(0.05, 0.1, 0.5))
)
tune_ctrl <- makeTuneControlGrid()
print(Sys.time())
tuned_result <- tuneParams(learner, task, resampling = makeResampleDesc("CV", iters = 5),
                           par.set = param_set, control = tune_ctrl)
print(Sys.time())
#svm_models <- list(); svm_rmses <- vector(); svm_rmsesIS <- vector()
i=3; print(Sys.time())
svm_models[[i]] <- svm(formula=as.formula(paste0("stage_train ~ lag_1+lag_2+lag_3+atmPres + gnossWind + 
    napaFlow + oceanWind + predNov")), 
                       data = dat_train_3, 
                       type= 'nu-regression',
                       kernel = "radial", # 'linear', 'sigmoid', 'polynomial', 'rbf' (radial)
                       nu = 0.5,
                       cost = 5, 
                       gamma = 0.1,
                       cross=0)
print(Sys.time())


set.seed(890)
print(Sys.time())
dat_train_0 <- dat_train %>% select(atmPres, gnossWind, napaFlow , oceanWind, predNov)
dat_train_0$stage <- stage_train
task <- makeRegrTask(data=dat_train_0, target='stage')
learner <- makeLearner('regr.svm', type='nu-regression')
param_set <- makeParamSet(
  makeDiscreteParam("kernel", values = c("radial")),#, "linear", "polynomial")),
  makeDiscreteParam("cost", values = c(1,2.5,5,10)),
  makeDiscreteParam("nu", values = c(0.5, 0.6)),
  makeDiscreteParam("gamma", values = c(0, 0.1))
  # makeDiscreteParam("kernel", values = c("radial")),#, "linear", "polynomial")),
  # makeDiscreteParam("cost", values = c(1, 2.5, 5)),
  # makeDiscreteParam("nu", values = c(0.1, 0.3, 0.5)),
  # makeDiscreteParam("gamma", values = c(0.05, 0.1, 0.5))
)
tune_ctrl <- makeTuneControlGrid()
print(Sys.time())
tuned_result <- tuneParams(learner, task, resampling = makeResampleDesc("CV", iters = 5),
                           par.set = param_set, control = tune_ctrl)
print(Sys.time())
#svm_models <- list(); svm_rmses <- vector(); svm_rmsesIS <- vector()
i=2; print(Sys.time())
svm_models[[i]] <- svm(formula=as.formula(paste0("stage_train~atmPres+gnossWind+napaFlow+oceanWind+predNov")), 
                       data = dat_train_0, 
                       type= 'nu-regression',
                       kernel = "radial", # 'linear', 'sigmoid', 'polynomial', 'rbf' (radial)
                       nu = 0.5,
                       cost = 10, 
                       gamma = 0.1,
                       cross=0)
print(Sys.time())

## DENSE FORECAST
modelNames <-  c('SVM24', 'SVM0', 'SVM3')
N_models=3;#length(svm_models); 
N_horizon=96; horizons <- c(1,24,48,72,96)
forecasts <-          lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1))); 
forecast_residuals <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1)));
forecasts_intervals <- lapply(1:N_models, function(x) list());
forecast_rmses <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon)));
forecast_maes <-  lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_bounds <-lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_horizon_rmses <- matrix(NA, nrow=N_models, ncol=length(horizons), dimnames = list(modelNames,  paste0(horizons,'hr'))); 
max_resids <- rep(NA, N_models);min_resids <- rep(NA, N_models);mean_rmses <- rep(NA, N_models);sd_rmses <- rep(NA, N_models);mean_maes <- rep(NA, N_models);mean_bounds <- rep(NA, N_models);
forecast_long_horizon <- list();forecast_residual_long_horizon <- list();forecasts_interval_long_horizon <- list();
forecast_rmse_long_horizon <- rep(NA,N_models);forecast_mae_long_horizon <- rep(NA,N_models);forecast_bound_long_horizon <- rep(NA,N_models);
N_lags <- c(24,0,3); 
start_tm <- Sys.time(); print(start_tm); set.seed(96)
## For each XGB model
m=1; # 

m=3
print(Sys.time())
dat_train_3 <- dat_train %>% select(lag_1, lag_2, lag_3, atmPres, gnossWind, napaFlow , oceanWind, predNov)
dat_train_3$stage <- stage_train
#svm_models <- list(); svm_rmses <- vector(); svm_rmsesIS <- vector()
i=3; print(Sys.time())
svm_models[[i]] <- svm(formula=as.formula(paste0("stage_train ~ lag_1+lag_2+lag_3+atmPres + gnossWind + 
    napaFlow + oceanWind + predNov")), 
                       data = dat_train_3, 
                       type= 'nu-regression',
                       kernel = "radial", # 'linear', 'sigmoid', 'polynomial', 'rbf' (radial)
                       nu = 0.5,
                       cost = 10, 
                       gamma = 0.1,
                       cross=0)
print(Sys.time())


test_dat <- as.matrix(as_tibble(dat_test) %>% dplyr::select(paste0('lag_',1:N_lags[m]), atmPres, gnossWind, napaFlow, oceanWind, predNov))
  ## Iterate through 2353 columns
for(fc in 1:(N_test-N_horizon+1)){
  print(paste0('fc: ',fc))
  ## set horizon window of 96 hours
  horizon <- fc:(fc+N_horizon-1)
  ## initialise lag responses at start of forecast window
  lag_window <- test_dat[fc,1:(N_lags[m])];
  ## iterate through 96 hours of forecast
  for(i in horizon){
    #print(paste0(m,'-',fc,'-',i))
    ## update lagged values for hour of interest + generate dmatrix of one row
    test_dat[i,1:(N_lags[m])] <- lag_window
    X_i <- matrix(test_dat[i,], ncol=dim(test_dat)[2], nrow=1, dimnames=list(i, colnames(test_dat))); 
    ## forecast hour of interest
    forecast_m <- predict(svm_models[[m]], X_i)
    forecasts[[m]][(i-fc+1),fc] <- forecast_m
    ## update lagged values with most recent forecast
    lag_window=c(forecast_m, lag_window)[1:(N_lags[m])];
  }
  ## calculate residuals and metrics for each 96 hour forecast
  resid_m <-  stage_test[horizon] - forecasts[[m]][,fc]
  forecast_residuals[[m]][,fc] <- resid_m
  forecast_rmses[[m]][fc] <- sqrt(mean(resid_m^2))
  forecast_maes[[m]][fc] <- mean(abs(resid_m))
  
}
## calculate aggregate metrics for model
for(h in 1:length(horizons)){
  forecast_horizon_rmses[m,h] <- sqrt(mean(forecast_residuals[[m]][horizons[h],]^2))
}
mean_rmses[m] <- round(mean(forecast_rmses[[m]]),6)
sd_rmses[m] <- round(sd(forecast_rmses[[m]]),6)
mean_maes[m] <- round(mean(forecast_maes[[m]]),6)
min_resids[m] <- round(min(forecast_residuals[[m]]),6)
max_resids[m] <- round(max(forecast_residuals[[m]]),6)
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

start_tm <- Sys.time(); print(start_tm); set.seed(96)
m=2
test_dat <- as.matrix(as_tibble(dat_test) %>% dplyr::select(atmPres, gnossWind, napaFlow, oceanWind, predNov))
## Iterate through 2353 columns
for(fc in 1:(N_test-N_horizon+1)){
  print(paste0(m,'-',fc))
  ## set horizon window of 96 hours
  horizon <- fc:(fc+N_horizon-1)
  ## iterate through 96 hours of forecast
  X_fc <- matrix(test_dat[horizon,], ncol=dim(test_dat)[2], nrow=N_horizon, dimnames=list(horizon, colnames(test_dat))); 
  ## forecast horizon
  forecast_m <- predict(svm_models[[m]], X_fc)
  forecasts[[m]][,fc] <- forecast_m
  ## calculate residuals and metrics for each 96 hour forecast
  resid_m <-  stage_test[horizon] - forecast_m
  forecast_residuals[[m]][,fc] <- resid_m
  forecast_rmses[[m]][fc] <- sqrt(mean(resid_m^2))
  forecast_maes[[m]][fc] <- mean(abs(resid_m))
  
}
## calculate aggregate metrics for model
for(h in 1:length(horizons)){
  forecast_horizon_rmses[m,h] <- sqrt(mean(forecast_residuals[[m]][horizons[h],]^2))
}
mean_rmses[m] <- round(mean(forecast_rmses[[m]]),6)
sd_rmses[m] <- round(sd(forecast_rmses[[m]]),6)
mean_maes[m] <- round(mean(forecast_maes[[m]]),6)
min_resids[m] <- round(min(forecast_residuals[[m]]),6)
max_resids[m] <- round(max(forecast_residuals[[m]]),6)
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))



highTideBin <- matrix(NA, nrow=96, ncol=2352)
for(t in 1:(2448-96)){
  highTideBin[,t] <- (stage_test[t:(t+95)]>0)
}
for(i in 1:length(forecast_residuals)){
  ht_rmse_sum <- c()
  for(j in 1:2352){ht_rmse_sum <- c(ht_rmse_sum, sqrt(mean((forecast_residuals[[i]][,j][highTideBin[,j]])^2)))}
  print(mean(ht_rmse_sum)*100)
  lt_rmse_sum <- c()
  for(j in 1:2352){lt_rmse_sum <- c(lt_rmse_sum, sqrt(mean((forecast_residuals[[i]][,j][!highTideBin[,j]])^2)))}
  print(mean(lt_rmse_sum)*100)
}



par(mfrow=c(1,2))
for(m in 1:2){
  xseq <- seq(from=-0.4, to=0.4, length=1000)
  yseq <- dnorm(xseq, mean=mean(forecast_residuals[[m]]), sd = sd(forecast_residuals[[m]]))
  yseq0 <- dnorm(xseq, mean=0, sd = sd(forecast_residuals[[m]]))
  hist(forecast_residuals[[m]],breaks=50, freq=F, main=modelNames[m], xlab='Forecast Residuals')
  lines(xseq, yseq, col='red')
  lines(xseq, yseq0, col='blue')
}



N_train <- dim(dat_train)[1]
for(m in 1:N_models){
  sigm <- sd(forecast_residuals[[m]])
  for(fc in 1:(N_test-N_horizon+1)){
    horizon=fc:(fc+95)
    t_value <- qt(0.975,df=N_train-1)
    pred <- forecasts[[m]][,fc]
    forecast_lower <- pred - (t_value+0.2) *sigm
    forecast_upper <- pred + (t_value+0.2) *sigm
    forecasts_intervals[[m]][[fc]] <- cbind(forecast_lower,forecast_upper)
    forecast_bounds[[m]][fc] <- round(length(intersect(which(stage_test[horizon]>forecast_lower),which(stage_test[horizon]<forecast_upper)))/N_horizon,4)
  }
  mean_bounds[m] <- round(mean(forecast_bounds[[m]]),6)
}
mean_bounds

forecast_hightide_rmses <- matrix(NA, ncol=N_models, nrow=N_horizon)
forecast_hightide_maes <-  matrix(NA, ncol=N_models, nrow=N_horizon)
greaterthan0 <- which(stage_test>0)
for(m in 1:N_models){
  for(h in 1:N_horizon){
    row_index <- h:(2352+h-1)
    print(range(greaterthan0+h-1))
    greater <-  tibble(gt0=greaterthan0-h+1) %>% filter(gt0>0) %>% filter(gt0<(N_test-N_horizon+1))
    print(range(greater$gt0))
    hightides <- forecast_residuals[[m]][h,greater$gt0]
    forecast_hightide_rmses[h,m] <- sqrt(mean(hightides^2))
    forecast_hightide_maes[h,m] <- mean(abs(hightides))
  }
}

m=2
horizon=1:N_test; test_dat <- as.matrix(as_tibble(dat_test) %>% dplyr::select(atmPres, gnossWind, napaFlow, oceanWind, predNov))
forecast_m <- predict(svm_models[[m]], test_dat)
forecast_long_horizon[[m]] <- forecast_m
t_value <- qt(0.975,df=N_train-1)
forecast_lower <- forecast_m - (t_value+0.2) *sigm
forecast_upper <- forecast_m + (t_value+0.2) *sigm
forecasts_interval_long_horizon[[m]] <- cbind(forecast_lower, forecast_upper)
resid_m <-  stage_test[horizon] - forecast_m
forecast_residual_long_horizon[[m]] <- resid_m
forecast_rmse_long_horizon[m] <- sqrt(mean(resid_m^2))
forecast_mae_long_horizon[m] <- mean(abs(resid_m))
forecast_bound_long_horizon[m] <- round(length(intersect(which(stage_test[horizon]>forecast_lower),which(stage_test[horizon]<forecast_upper)))/N_test,4)


t(t(colMeans(forecast_hightide_rmses)))
t(t(colMeans(forecast_hightide_maes)))
t(forecast_hightide_rmses[horizons,])

sqrt(mean((svm_models[[2]]$residuals)^2))

forecast_horizon_rmses*100
#cbind(mean_rmses, sd_rmses, mean_maes)
t(t(mean_rmses*100))
t(t(sd_rmses*100))
t(t(mean_maes*100))
#t(t(mean_bounds))
t(t(min_resids*100))
t(t(max_resids*100))