rm(list=ls())


# Data Directory
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")

setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='xgboost_stage_dense2.RData')
# save.image(file='xgboost_stage_dense_importance.RData')
#load(file = 'xgboost_stage.RData')
#load(file = 'xgboost_stage_ext.RData')
#load(file='xgboost_stage_dense.RData')
#load(file='xgboost_stage_dense2.RData')


#Libraries
#library(astsa)
library(forecast)
library(ggplot2)
library(gridExtra)
library(nlme)
library(scales)
library(tidyverse)
library(tseries)
#library(urca)
#library(vars)
library(xgboost)
#library(zoo)

# Load Data
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)
# Standardise - centre and scale data
#X[,2:5] <- scale(X[,2:5], center = T, scale=T)
#atmP <- (X$AtmPres-min(X$AtmPres))/max(X$AtmPres-min(X$AtmPres)); 
for(i in 2:5){X[,i] <- rescale(X[,i], to=c(-1,1))}
# Format datetime
X$datetime <- as.POSIXct(X$datetime, tz="GMT")
y_nov$datetime <- as.POSIXct(y_nov$datetime, tz="GMT")

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

#a <- acf(dat_nov$stageNov, lag.max = 720)
#p <- pacf(dat_nov$stageNov, lag.max=96)

# Create time series object
N <- dim(X)[1]; N_test <- 2448; N_train <- N-N_test; N_lags=24
dat_lags <- matrix(NA, nrow=N, ncol=N_lags)
for(i in 1:N_lags){dat_lags[(1+i):(N),i]=dat_nov$stageNov[(1):(N-i)]}
colnames(dat_lags) <- paste0('lag_',1:N_lags)

dat_train <- na.omit(cbind(dat_lags[1:N_train,], dat_nov[1:N_train,]))
dat_test <- cbind(dat_lags[(1+N_train):N,], dat_nov[(1+N_train):N,])

dt_train <- dat_nov$datetime[(1+N_lags):N_train]
dt_test <- dat_nov$datetime[(1+N_train):N]


stage_train <- dat_train$stageNov
stage_test <- dat_test$stageNov
dat_train <- na.omit(ts(data=(dat_train%>%select(colnames(dat_lags), atmPres, gnossWind, napaFlow, oceanWind, predNov)), start=c(2019, 20+N_lags),
                end=c(2022,(6477.562-N_test)), frequency=365.24219*24, class="matrix"))
dat_test <- na.omit(ts(data=(dat_test%>%select(colnames(dat_lags), atmPres, gnossWind, napaFlow, oceanWind, predNov)), start=c(2022,(6477.562-N_test+1)),
               end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix"))

dtrain <- list(xgb.DMatrix(data = as.matrix(dat_train), label = stage_train),
               xgb.DMatrix(data = as.matrix(dat_train[,c(1:3, 25:29)]), label = stage_train),
               xgb.DMatrix(data = as.matrix(dat_train[,c(25:29)]), label = stage_train))
dtest <- list(xgb.DMatrix(data = as.matrix(dat_test)),# label = resNov_test),
              xgb.DMatrix(data = as.matrix(dat_test[,c(1:3, 25:29)])),# label = resNov_test),
              xgb.DMatrix(data = as.matrix(dat_test[,c(25:29)])))#, label = resNov_test))



dtrain <- list(xgb.DMatrix(data = as.matrix(dat_train), label = stage_train),
               xgb.DMatrix(data = as.matrix(dat_train[,c(1:3, 25:29)]), label = stage_train),
               xgb.DMatrix(data = as.matrix(dat_train[,c(25:29)]), label = stage_train))
dtest <- list(xgb.DMatrix(data = as.matrix(dat_test)),# label = resNov_test),
              xgb.DMatrix(data = as.matrix(dat_test[,c(1:3, 25:29)])),# label = resNov_test),
              xgb.DMatrix(data = as.matrix(dat_test[,c(25:29)])))#, label = resNov_test))

 
param_grid <- expand.grid(
  eta = c(0.1),#, 0.5, 1),                   # step size / learning rate
  max_depth = c(6, 9, 12),                    # Larger -> greater complexity + risk of overfit
  subsample = c(0.6, 0.8, 1),              # Subsample ratio (proportiton of training data used) prevents overfitting
  colsample_bytree = c(0.8, 1),  #0.6,       # subsample ratio of columns when constructing each tree (one subsample per tree)
  #min_child_weight = c(1, 3, 5),             # Larger -> More conservative (Less partitioning)
  gamma = c(0, 0.1),#, 0.5),                    # min loss reduction to make partition: Larger -> More conservative (Less partitioning)
  lambda=c(0, 1, 5)                           # regularisation penalty
)

cv_results <- list(lag24=data.frame(), lag3=data.frame(), lag0=data.frame())

set.seed(4321)
for(m in 1:2){
  print(paste0('Model: ',m,' - ', Sys.time()))
for(i in 1:nrow(param_grid)) {
  print(paste0('Iteration: ',i,' - ', Sys.time()))
  params <- list(
    objective = "reg:squarederror",
    booster = "gbtree",
    eta = param_grid$eta[i],
    max_depth = param_grid$max_depth[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma = param_grid$gamma[i],
    lambda=param_grid$lambda[i]
  )
  
  xgb_cv <- xgb.cv(
    params = params,
    data = dtrain[[m]],
    nrounds = 200,
    nfold = 10,
    early_stopping_rounds = 50,
    verbose = 0
  )
  
  cv_results[[m]] <- rbind(cv_results[[m]], data.frame(
    eta = param_grid$eta[i],
    max_depth = param_grid$max_depth[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma = param_grid$gamma[i],
    lambda=param_grid$lambda[i],
    best_nrounds = xgb_cv$best_iteration,
    best_rmse_IS = min(xgb_cv$evaluation_log$train_rmse_mean),
    best_rmse_OOS = min(xgb_cv$evaluation_log$test_rmse_mean)
  ))
  
  saveRDS(cv_results, file = "cv_results_intermediate.rds")
}
}

# cv_results <- readRDS("cv_results_intermediate.rds")
#loaded_data <- readRDS("filename.rds")

cv_results[[3]] %>% arrange(best_rmse_OOS)

#write.csv(cv_results[[1]], file='xgboost_stage_params.csv', col.names=T)

param_grid2 <- expand.grid(
  eta = c(0.1, 0.2),                   # step size / learning rate
  max_depth = c(9),                    # Larger -> greater complexity + risk of overfit
  subsample = c(0.8, 0.9, 1),              # Subsample ratio (proportiton of training data used) prevents overfitting
  colsample_bytree = c(0.8),       # subsample ratio of columns when constructing each tree (one subsample per tree)
  #min_child_weight = c(1, 3, 5),             # Larger -> More conservative (Less partitioning)
  gamma = c(0),                    # min loss reduction to make partition: Larger -> More conservative (Less partitioning)
  lambda=c(0,2,5)
)


param_grid2 <- expand.grid(
  eta = c(0.1, 0.2),                   # step size / learning rate
  max_depth = c(9),                    # Larger -> greater complexity + risk of overfit
  subsample = c(0.8, 0.9, 1),              # Subsample ratio (proportiton of training data used) prevents overfitting
  colsample_bytree = c(0.8),       # subsample ratio of columns when constructing each tree (one subsample per tree)
  #min_child_weight = c(1, 3, 5),             # Larger -> More conservative (Less partitioning)
  gamma = c(0),                    # min loss reduction to make partition: Larger -> More conservative (Less partitioning)
  lambda=c(0,2,5)
)


cv_results_detail <- list(data.frame(),data.frame())

set.seed(4321); m=2
for (i in 1:nrow(param_grid2)) {
  print(paste0('Iteration: ',i,' - ', Sys.time()))
  params <- list(
    objective = "reg:squarederror",
    booster = "gbtree",
    eta = param_grid2$eta[i],
    max_depth = param_grid2$max_depth[i],
    subsample = param_grid2$subsample[i],
    colsample_bytree = param_grid2$colsample_bytree[i],
    gamma = param_grid2$gamma[i],
    lambda=param_grid2$lambda[i]
  )
  
  xgb_cv <- xgb.cv(
    params = params,
    data = dtrain[[m]],
    nrounds = 200,
    nfold = 10,
    early_stopping_rounds = 50,
    verbose = 0
  )
  
  cv_results_detail[[m]] <- rbind(cv_results_detail[[m]], data.frame(
    eta = param_grid2$eta[i],
    max_depth = param_grid2$max_depth[i],
    subsample = param_grid2$subsample[i],
    colsample_bytree = param_grid2$colsample_bytree[i],
    gamma = param_grid2$gamma[i],
    lambda=param_grid2$lambda[i],
    best_nrounds = xgb_cv$best_iteration,
    best_rmse_IS = min(xgb_cv$evaluation_log$train_rmse_mean),
    best_rmse_OOS = min(xgb_cv$evaluation_log$test_rmse_mean)
  ))
  
  #saveRDS(cv_results_detail, file = "cv_results_detail_intermediate.rds")
}

cv_results_detail[[1]] %>% arrange(best_rmse_OOS)
cv_results_detail[[2]] %>% arrange(best_rmse_OOS)

params24 <- list(
  objective = "reg:squarederror",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 9,
  subsample = 0.8,
  colsample_bytree =0.8,
  gamma = 0,
  lambda=5
)
params3 <- list(
  objective = "reg:squarederror",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 9,
  subsample = 0.8,
  colsample_bytree =0.8,
  gamma = 0,
  lambda=5
)

params0 <- list(
  objective = "reg:squarederror",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 9,
  subsample = 1.0,
  colsample_bytree =1.0,
  gamma = 0,
  lambda=5
)

xgb_cv_24 <- xgb.cv(
  params = params24,
  data = dtrain[[1]],
  nrounds = 500,
  nfold = 10,
  early_stopping_rounds = 50,
  verbose = 0
)
xgb_cv_3 <- xgb.cv(
  params = params3,
  data = dtrain[[2]],
  nrounds = 500,
  nfold = 10,
  early_stopping_rounds = 50,
  verbose = 0
)

xgb_models <- list(); xgb_forecasts <- list(); params <- list(params24, params3); 
test_rmses <- rep(NA, 2); test_maes <- rep(NA, 2); xgb_resids <- list(); max_resids <- rep(NA,2)

set.seed(4321)
for(i in 3){
  xgb_models[[i]] <- xgb.train(
    params = params0,
    data = dtrain[[i]],
    nrounds = 1000,
    feval = NULL,
    eval_metric = "rmse"
  )
  #xgb_models_lower[[i]] <- xgb.train
}


IS_pred <- list(); IS_resid <- list(); IS_RMSE <- vector(length=N_models)
for(m in 1:3){
  IS_pred[[m]] <- predict(xgb_models[[m]], dtrain[[m]])
  IS_resid[[m]] <-  stage_train - IS_pred[[m]]
  IS_RMSE[m] <- sqrt(mean(IS_resid[[m]]^2))
}

forecasts[[3]] <- matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1))
forecast_residuals[[3]] <- matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1))
forecast_rmses[[3]] <- rep(NA, (N_test - N_horizon))
forecast_maes[[3]]<- rep(NA, (N_test - N_horizon))
forecast_bounds[[3]]<- rep(NA, (N_test - N_horizon))
forecast_horizon_rmses <- rbind(forecast_horizon_rmses, rep(NA, 5)); rownames(forecast_horizon_rmses) <- modelNames
max_resids <- c(max_resids, NA)
min_resids <- c(min_resids, NA)
mean_rmses <- c(mean_rmses, NA)
sd_rmses <- c(sd_rmses, NA)
mean_maes <- c(mean_maes, NA)
mean_bounds <- c(mean_bounds, NA)
forecast_hightide_rmses <- cbind(forecast_hightide_rmses, rep(NA, N_horizon))
forecast_hightide_maes <- cbind(forecast_hightide_maes, rep(NA, N_horizon))

modelNames <-  c('XGB24', 'XGB3', 'XGB0')
N_models=length(xgb_models); N_horizon=96; horizons <- c(1,24,48,72,96)
forecasts <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1))); 
forecast_residuals <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1)));
forecasts_intervals <- lapply(1:N_models, function(x) list());
forecast_rmses <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon)));
forecast_maes <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_bounds <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_horizon_rmses <- matrix(NA, nrow=N_models, ncol=5, dimnames = list(modelNames,  paste0(horizons,'hr'))); 
max_resids <- rep(NA, N_models);min_resids <- rep(NA, N_models);mean_rmses <- rep(NA, N_models);sd_rmses <- rep(NA, N_models);mean_maes <- rep(NA, N_models);mean_bounds <- rep(NA, N_models);
forecast_long_horizon <- list();forecast_residual_long_horizon <- list();forecasts_interval_long_horizon <- list();
forecast_rmse_long_horizon <- rep(NA,N_models);forecast_mae_long_horizon <- rep(NA,N_models);forecast_bound_long_horizon <- rep(NA,N_models);
N_lags <- c(24,3,0); 
start_tm <- Sys.time(); print(start_tm); set.seed(96)
## For each XGB model
for(m in 3){#2:N_models){
  test_dat <- as.matrix(as_tibble(dat_test) %>% dplyr::select(paste0('lag_',1:N_lags[m]), atmPres, gnossWind, napaFlow, oceanWind, predNov))
  ## Iterate through 2353 columns
  for(fc in 1:(N_test-N_horizon+1)){
    ## set horizon window of 96 hours
    horizon <- fc:(fc+N_horizon-1)
    ## initialise lag responses at start of forecast window
    lag_window <- test_dat[fc,1:(N_lags[m])];
    ## iterate through 96 hours of forecast
    for(i in horizon){
      print(paste0(m,'-',fc,'-',i))
      ## update lagged values for hour of interest + generate dmatrix of one row
      test_dat[i,1:(N_lags[m])] <- lag_window
      X_i <- matrix(test_dat[i,], ncol=dim(test_dat)[2], nrow=1, dimnames=list(i, colnames(test_dat))); 
      dX_i <- xgb.DMatrix(data = X_i)
      ## forecast hour of interest
      forecast_m <- predict(xgb_models[[m]], dX_i)
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
  
}
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

par(mfrow=c(1,2))
for(m in 1:2){
xseq <- seq(from=-0.4, to=0.4, length=1000)
yseq <- dnorm(xseq, mean=mean(forecast_residuals[[m]]), sd = sd(forecast_residuals[[m]]))
yseq0 <- dnorm(xseq, mean=0, sd = sd(forecast_residuals[[m]]))
hist(forecast_residuals[[m]],breaks=50, freq=F, main=modelNames[m], xlab='Forecast Residuals')
lines(xseq, yseq, col='red')
lines(xseq, yseq0, col='blue')
}

round(sqrt(rowMeans(forecast_residuals[[2]]^2)[1:96])*100,1)

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

t(t(colMeans(forecast_hightide_rmses)))
t(t(colMeans(forecast_hightide_maes)))
t(forecast_hightide_rmses[horizons,])




par(mfrow=c(1,1))

m=2
fc=647
plotint=fc:(fc+95)
plot(stage_test[plotint]~dt_test[plotint], type='l',lwd=3,col='grey', ylim=range(unlist(forecasts)))
lines(forecasts[[m]][,fc]~dt_test[plotint], col='blue',lwd=2)
lines(forecasts_intervals[[m]][[fc]][,1]~dt_test[plotint], col='red')
lines(forecasts_intervals[[m]][[fc]][,2]~dt_test[plotint], col='red')

m=1
fc=647
plotint=fc:(fc+95)
plot(stage_test[plotint]~dt_test[plotint], type='l',lwd=3,col='grey', ylim=range(unlist(forecasts)))
lines(forecasts[[m]][,fc]~dt_test[plotint], col='blue',lwd=2)
lines(forecast_lower~dt_test[plotint], col='red')
lines(forecast_upper~dt_test[plotint], col='red')


forecast_horizon_rmses
t(t(mean_rmses))
t(t(sd_rmses))
t(t(mean_maes))
t(t(min_resids))
t(t(max_resids))

[1,] 0.074745
[2,] 0.085068
[3,] 0.092077

predictions <- xgb_forecasts[[1]]
xgb_resid <- stage_test-predictions
rmse <- sqrt(mean((xgb_resid)^2))
mae <- mean(abs(xgb_resid))
sigm <- sd(xgb_resid)
N_train <- dim(dat_train)[1]
t_value <- qt(0.975,df=N_train-1)

pred_lower <- predictions - t_value * sqrt(2*sigm^2)
pred_upper <- predictions + t_value * sqrt(2*sigm^2)

propBound <- mean(ifelse((stage_test>=pred_lower),1,0)*ifelse((stage_test<=pred_upper),1,0))

max(xgb_resid)
min(xgb_resid)

i=1;
model_fit <- predict(xgb_models[[i]], dtrain[[i]])
sqrt(mean((stage_train-model_fit)^2))

# set.seed(4321)
# n_models <- 30
# boot_preds <- matrix(NA, nrow = nrow(dat_test), ncol = n_models)
# i=1
# for (b in 1:n_models) {
#   print(paste0('Bootstrap Iter: ',b))
#   boot_sample <- sample(1:(N_train-N_lags), replace = TRUE)
#   dboot <- xgb.DMatrix(data = as.matrix(dat_train[boot_sample,]), label = stage_train[boot_sample])
#   
#   params <- list(
#     objective = "reg:squarederror",
#     booster = "gbtree",
#     eta = 0.1,
#     max_depth = 9,
#     subsample = 0.8,
#     colsample_bytree =0.8,
#     gamma = 0,
#     lambda=5
#   )
#   
#   m_boot <- xgb.train(
#     params = params,
#     data = dboot,
#     nrounds = 100,
#     feval = NULL,
#     eval_metric = "rmse"
#   )
#   
#   boot_preds[, b] <- predict(m_boot, dtest[[i]])
# }
# 
# lower_bound <- apply(boot_preds, 1, quantile, probs = 0.025)
# upper_bound <- apply(boot_preds, 1, quantile, probs = 0.975)
# mean_prediction <- apply(boot_preds, 1, mean)
# 
# # Combine the results into a data frame
# results <- data.frame(
#   Mean = mean_prediction,
#   Lower = lower_bound,
#   Upper = upper_bound
# )

plot_n=72
plot(stage_test~dt_test, type='l', lwd=4, col='darkgrey', 
     main=paste0('Test Forecast: xgb',i, ' - Full Test Range'), ylab='stageNov', xlab='date')
lines(xgb_forecasts[[1]]~dt_test, col='black', lty=1, lwd=2)
# lines(results$Mean~dt_test, col=1, lty=1, lwd=2)
# lines(results$Lower~dt_test, col=2, lty=1, lwd=2)
# lines(results$Upper~dt_test, col=2, lty=1, lwd=2)

plot(stage_test[1:plot_n]~dt_test[1:plot_n], type='l', lwd=4, col='darkgrey', 
     main=paste0('Test Forecast: xgb',i, ' - ', plot_n,' hours'), ylab='stageNov', xlab='date')
lines(predictions[1:plot_n]~dt_test[1:plot_n], col=1, lty=1, lwd=2)
lines(pred_lower[1:plot_n]~dt_test[1:plot_n], col=2, lty=1, lwd=2)
lines(pred_upper[1:plot_n]~dt_test[1:plot_n], col=2, lty=1, lwd=2)


set.seed(4321)
for(i in 1:2){
  xgb_models[[i]] <- xgb.train(
    params = params[[i]],
    data = dtrain[[i]],
    nrounds = 1000,                    # Number of trees
    watchlist = list(train = dtrain[[i]]),
    print_every_n = 50,
    early_stopping_rounds = 50        # Stop if no improvement
  )
  
  xgb_forecasts[[i]] <- predict(xgb_models[[i]], newdata=dtest[[i]])
  xgb_resids[[i]] <- stage_test - xgb_forecasts[[i]]
}

par(mfrow=c(2,4))
for(i in 1:length(xgb_models)){
  plot(stage_test~dt_test, type='l', lwd=4, col='darkgrey', 
       main=paste0('Test Forecast: xgb',i, ' - Full Test Range'), ylab='stageNov', xlab='date')
  lines(xgb_forecasts[[i]]~dt_test, col=1, lty=1, lwd=2)
  
  plot(xgb_resids[[i]]~dt_test, type='l', lwd=2, ylab='modelResidual', xlab='date')
  
  test_rmses[i] <- sqrt(mean((xgb_resids[[i]])^2))
  test_maes[i] <- mean(abs(xgb_resids[[i]]))
  max_resids[i] <- max(xgb_resids[[i]])
  
  plot_n=72
  plot(stage_test[1:plot_n]~dt_test[1:plot_n], type='l', lwd=4, col='darkgrey', 
       main=paste0('Test Forecast: xgb',i, ' - ', plot_n,' hours'), ylab='stageNov', xlab='date')
  lines(xgb_forecasts[[i]][1:plot_n]~dt_test[1:plot_n], col=1, lty=1, lwd=2)
  
  plot(xgb_resids[[i]][1:plot_n]~dt_test[1:plot_n], type='l', lwd=2, ylab='modelResidual', xlab='date')
}
print(cbind(test_rmses, test_maes, max_resids))








start_tm <- Sys.time(); print(start_tm); set.seed(96)
for(m in 3){#2:N_models){
  test_dat <- as.matrix(as_tibble(dat_test) %>% dplyr::select(atmPres, gnossWind, napaFlow, oceanWind, predNov))
  ## Iterate through 2353 columns
  for(fc in 1:(N_test-N_horizon+1)){
    ## set horizon window of 96 hours
    horizon <- fc:(fc+N_horizon-1)
    X_fc <- matrix(test_dat[horizon,], ncol=dim(test_dat)[2], nrow=N_horizon, dimnames=list(horizon, colnames(test_dat))); 
    dX_fc <- xgb.DMatrix(data = X_fc)
    ## forecast hour of interest
    forecast_m <- predict(xgb_models[[m]], dX_fc)
    forecasts[[m]][,fc] <- forecast_m
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
  
  sigm <- sd(forecast_residuals[[m]])
  for(fc in 1:(N_test-N_horizon+1)){
    horizon=fc:(fc+95)
    t_value <- qt(0.975,df=N_train-1)
    pred <- forecasts[[m]][,fc]
    forecast_lower <- pred - (t_value+0.2) *sigm
    forecast_upper <- pred + (t_value+0.2) *sigm
    forecasts_intervals0[[fc]] <- cbind(forecast_lower,forecast_upper)
    forecast_bounds[[m]][fc] <- round(length(intersect(which(stage_test[horizon]>forecast_lower),which(stage_test[horizon]<forecast_upper)))/N_horizon,4)
  }
  mean_bounds[m] <- round(mean(forecast_bounds[[m]]),6)
  
  
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
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

[1,] 0.07757898 0.07608826 0.07589380 0.07564546 0.07579265
[2,] 0.08974453 0.08946248 0.08962108 0.08963120 0.08966090
[3,] 0.08736086 0.08716535 0.08739530 0.08736437 0.08762289


rmse <- function(x){
  return(sqrt(mean((x)^2)))
}

plot(apply(forecast_residuals[[1]][,800:1200], 1, rmse))
plot(forecast_rmses[[1]])
mean(apply(forecast_residuals[[1]], 1, rmse))
mean(forecast_rmses[[1]])
library(pdp)
########
### Partial Dependence Plots
########

pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/96-hour Forecasts/',modelNames[i],'-',n_forecast,'.pdf'), width = 9, height = 4)
par(mfrow=c(1,1), mar=c(4,4,2,0.5)); 

dev.off()

# atmPres
pdp.AP.24 <- partial(object = xgb_models[[1]], pred.var = "atmPres", train = dat_train)
names(pdp.AP.24)[1] <- "x"
pdp.AP.24$model <- "XGB24" 
pdp.AP.3 <- partial(object = xgb_models[[2]], pred.var = "atmPres", train = dat_train[,c(1:3,25:29)])
names(pdp.AP.3)[1] <- "x"
pdp.AP.3$model <- "XGB3"
pdp.AP.0 <- partial(object = xgb_models[[3]], pred.var = "atmPres", train = dat_train[,25:29])
names(pdp.AP.0)[1] <- "x"
pdp.AP.0$model <- "XGB0"
pdp.AP.all <- rbind(pdp.AP.24, pdp.AP.3, pdp.AP.0)

# Tidal Pred
pdp.TP.24 <- partial(object = xgb_models[[1]], pred.var = "predNov", train = dat_train)
names(pdp.TP.24)[1] <- "x"
pdp.TP.24$model <- "XGB24" 
pdp.TP.3 <- partial(object = xgb_models[[2]], pred.var = "predNov", train = dat_train[,c(1:3,25:29)])
names(pdp.TP.3)[1] <- "x"
pdp.TP.3$model <- "XGB3"
pdp.TP.0 <- partial(object = xgb_models[[3]], pred.var = "predNov", train = dat_train[,25:29])
names(pdp.TP.0)[1] <- "x"
pdp.TP.0$model <- "XGB0"
pdp.TP.all <- rbind(pdp.TP.24, pdp.TP.3, pdp.TP.0)

# Lag
pdp.l1.24 <- partial(object = xgb_models[[1]], pred.var = "lag_1", train = dat_train)
names(pdp.l1.24)[1] <- "x"
pdp.l1.24$model <- "XGB24" 
pdp.l1.3 <- partial(object = xgb_models[[2]], pred.var = "lag_1", train = dat_train[,c(1:3,25:29)])
names(pdp.l1.3)[1] <- "x"
pdp.l1.3$model <- "XGB3"
pdp.l1.0 <- pdp.l1.3[1,]; pdp.l1.0$yhat <- 10; pdp.l1.0$model <- 'XGB0'

pdp.l.all <- rbind(pdp.l1.24, pdp.l1.3)

yrange <- range(c(pdp.AP.all$yhat,pdp.TP.all$yhat,pdp.l.all$yhat))

pdp.l.all <- rbind(pdp.l.all,pdp.l1.0)

pdp.AP <- ggplot(pdp.AP.all, aes(x = x, y = yhat, color = model)) +
  geom_line(linewidth = 1) +  # use 'linewidth' instead of 'size'
  scale_color_manual(values = c("XGB24" = "red", 
                                "XGB3"  = "blue", 
                                "XGB0"  = "green")) +
  labs(x = "AP (scaled [-1,1])",
       y = "Partial Dependence of Stage (cm)",
       title = "",
       color = "Model") +
  ylim(yrange) +
  theme_minimal() +
  theme(
    #panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # add border
    legend.position = "none"
  )

pdp.TP <- ggplot(pdp.TP.all, aes(x = x, y = yhat, color = model)) +
  geom_line(linewidth = 1) +  # use 'linewidth' instead of 'size'
  scale_color_manual(values = c("XGB24" = "red", 
                                "XGB3"  = "blue", 
                                "XGB0"  = "green")) +
  labs(x = "TP",
       y = "",
       title = "",
       color = "Model") +
  ylim(yrange) +
  theme_minimal() +
  theme(
    #panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # add border
    legend.position = "none"
  )


pdp.lag <- ggplot(pdp.l.all, aes(x = x, y = yhat, color = model)) +
  geom_line(linewidth = 1) +  # use 'linewidth' instead of 'size'
  scale_color_manual(values = c("XGB24" = "red", 
                                "XGB3"  = "blue", 
                                "XGB0"  = "green")) +
  labs(x = "Lag 1",
       y = "",
       title = "",
       color = "Model") +
  ylim(yrange) +
  theme_minimal() +
  theme(
    #panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # add border
  )

pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/PDP/xgb_pdp.pdf'), 
    width = 9, height = 3)
#par(mfrow=c(1,1), mar=c(4,4,2,0.5)); 
grid.arrange(pdp.AP, pdp.TP, pdp.lag, nrow = 1, widths = c(1, 1, 1.38))
dev.off()



highTideBin <- matrix(NA, nrow=96, ncol=2352)
for(t in 1:(2448-96)){
  highTideBin[,t] <- (stage_test[t:(t+95)]>0)
}
for(i in 1:length(forecast_residuals)){
  ht_rmse_sum <- c()
  for(j in 1:2352){ht_rmse_sum <- c(ht_rmse_sum, sqrt(mean((forecast_residuals[[i]][,j][highTideBin[,j]])^2)))}
  print(mean(ht_rmse_sum)*100)
  lt_rmse_sum <- c()
  for(j in 1:2352){lt_rmse_sum <- c(lt_rmse_sum,sqrt(mean((forecast_residuals[[i]][,j][!highTideBin[,j]])^2)))}
  print(mean(lt_rmse_sum)*100)
}

error_margins <- c(0.05,0.1,0.15,0.2, 0.25)
for(m in 1:length(forecast_residuals)){
  
  for(c in 1:5){
    print(paste0(modelNames[m],'-Tot-',error_margins[c],': ',100*round(sum(abs(forecast_residuals[[m]])<=error_margins[c])/prod(dim(forecast_residuals[[m]])),3)))
    #pos_stage_count <- 0
    #for(j in 1:2352){pos_stage_count <- pos_stage_count+sum(abs(forecast_residuals[[m]][,j][highTideBin[,j]])<=error_margins[c])}
    #print(paste0(modelNames[m],'-Pos-',error_margins[c],': ',100*round(pos_stage_count/sum(colSums(highTideBin)),3)))
  }
  
  #lt_rmse_sum <- 0
  #for(j in 1:2352){lt_rmse_sum <- lt_rmse_sum+sqrt(mean((forecast_residuals[[i]][,j][!highTideBin[,j]])^2))}
  #print(lt_rmse_sum/2352*100)
}

m=1
round(sum(abs(forecast_residuals[[m]])<=0.05)/prod(dim(forecast_residuals[[m]])),3)
round(sum(abs(forecast_residuals[[m]])<=0.1)/prod(dim(forecast_residuals[[m]])),3)
round(sum(abs(forecast_residuals[[m]])<=0.15)/prod(dim(forecast_residuals[[m]])),3)
round(sum(abs(forecast_residuals[[m]])<=0.2)/prod(dim(forecast_residuals[[m]])),3)

(lt_rmse_sum + ht_rmse_sum)/2352*100/2

ht_rmse_sum <- 0
for(j in 1:2352){ht_rmse_sum <- ht_rmse_sum+sqrt(mean((forecast_residuals[[i]][,j])^2))}
print(ht_rmse_sum/2352*100)

print(sqrt(sum((forecast_residuals[[i]])^2)/96)/2353*100)

n_forecast=914; horizon=(n_forecast):(n_forecast+95); i=1#2205 #420
rng=range(c(#forecasts[[1]][,n_forecast],
  #forecasts[[2]][,n_forecast],
  forecasts[[i]][,n_forecast],stage_test[horizon]))
#forecasts[[4]][,n_forecast],
#forecasts[[5]][,n_forecast],
#forecasts[[6]][,n_forecast]),

pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/96-hour Forecasts/',modelNames[i],'-',n_forecast,'.pdf'), width = 9, height = 4)
par(mfrow=c(1,1), mar=c(4,4,2,0.5)); 
plot(stage_test[horizon]~dt_test[horizon], col='darkgrey', pch=19, cex=1.2,
     main=paste0(modelNames[i],' 96-hr Forecast starting ',dt_test[horizon][1]), 
     ylab='Stage (m)',xlab='Datetime',
     ylim=rng)
lines(forecasts[[i]][,n_forecast]~dt_test[horizon], col='blue', lwd=2)
dev.off()

xgb_importance <- list();  par(mfrow=c(3,1))
for(i in 1:3){
  xgb_importance[[i]] <- xgb.importance(model=xgb_models[[i]])
  xgb_importance[[i]][,2:4] <- round(xgb_importance[[i]][,2:4],3)
  print(xgb_importance[[i]]) #%>% arrange(desc(Cover)))
  xgb.plot.importance(xgb_importance[[i]], main=paste0('XGB',N_lags[i]))
}

xgb_importance[[1]] %>% as.tibble() %>% left_join(as.tibble(xgb_importance[[2]]), join_by(Feature==Feature)) %>% 
  left_join(as.tibble(xgb_importance[[3]]), join_by(Feature==Feature)) %>% 
  mutate(Feature = factor(Feature, levels = c('predNov','atmPres','gnossWind','napaFlow','oceanWind', lag_names))) %>% 
  arrange(Feature) %>% select(Feature, Importance_24=Importance.x, Importance_3=Importance.y, Importance_NoLag=Importance,
                              Cover_24=Cover.x, Cover_3=Cover.y, Cover_NoLag=Cover,
                              Frequency_24=Frequency.x, Frequency_3=Frequency.y, Frequency_NoLag=Frequency) %>%
  #arrange(desc(Importance_3)) %>%
  print(n=29) 


# Gain: Measures the contribution of each feature to reducing the loss function across all splits in the trees.
#       Greater value indicates more critical role in model performance
# Cover: Indicates the relative proportion of observations (or "coverage") that a feature contributes to in the splits.
#        Greater value indicates feature split more observations
# Frequency: Indicates relative frequency of appearance of the feature in the decision trees.
#            Greater value indicates greater presence/increased use, even if splits do not contribute to reducing error.




