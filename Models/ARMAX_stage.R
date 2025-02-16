rm(list=ls())


# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='ARMIAX_stage.RData') # save.image(file='ARMAX_train.RData')  save.image(file='ARMAX_seas_temp.RData')
#load(file='ARMIAX_stage.RData')
# save.image(file='ARMIAX_stage_denseforecast.RData') 
# load(file='ARMIAX_stage_denseforecast.RData') 



#Libraries
library(astsa)
library(forecast)
library(ggplot2)
library(gridExtra)
library(nlme)
library(scales)
library(tidyverse)
library(tseries)
library(urca)
library(vars)
library(zoo)

# Load Data
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)
# Standardise - centre and scale data
#X[,2:5] <- scale(X[,2:5], center = T, scale=T)
#atmP <- (X$AtmPres-min(X$AtmPres))/max(X$AtmPres-min(X$AtmPres)); 
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

# include seasonal terms for arimax; use stepwise parameters selection using AIC 
# Create time series object
#dt <- dat_nov[,1]
#dat <- ts(data=dat_nov[2:6], start=c(2019, 20),
#          end=c(2022,6477.562), frequency=365.24219*24, class="matrix") #end=c(2021,12,28,15)
#tm_yr <- time(dat); start_tm <- min(tm_yr); end_tm <- max(tm_yr); yrs <- end_tm-start_tm; tm_yr <- ((tm_yr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#tm_28 <- time(dat); start_tm <- min(tm_28); end_tm <- max(tm_28); yrs <- (end_tm-start_tm)*365.24219/27.33; tm_28 <- ((tm_28-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#dat <- as.data.frame(cbind(sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:4] <- c('sin_28','cos_28','sin_yr','cos_yr')


# Create time series object
N <- dim(X)[1]; N_test <- 2448; N_train <- N-N_test
dat <- dat_nov %>% dplyr::select(datetime, atmPres, gnossWind, napaFlow, oceanWind, predNov, stageNov)
dt_train <- dat$datetime[1:N_train]
dt_test <- dat$datetime[(1+N_train):N]
dat_train <- ts(data=dat[1:N_train,2:7], start=c(2019, 20),
          end=c(2022,(6477.562-N_test)), frequency=365.24219*24, class="matrix") 
dat_test <- as.data.frame(ts(data=dat[(1+N_train):N,2:7], start=c(2022,(6477.562-N_test+1)),
               end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix")) 
stage_test <- dat_test$stageNov; dat_test <- as.data.frame(dat_test[,-6]);
dat_train <- as.data.frame(dat_train);

kpss.test(dat_train$stageNov)                     # Tests for stationarity
kpss.test(diff(dat_train$stageNov))                     # Tests for stationarity
Box.test(dat_train$stageNov, type = "Ljung-Box")  # Tests for autocorrelation

# pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/ACF.pdf'), 
#     width = 9, height = 3)
# par(mfrow=c(1,2), mar=c(4,4,2,1))
# pa <- pacf(dat_train$stageNov, lag.max = 96, main=''); #abline(h=c(-0.05, 0.05), lty=2, col='green');abline(h=c(-0.1, 0.1), lty=2, col='red')
# title(main='Stage PACF (96 hours)', line = 0.5)
# a <- acf(dat_train$stageNov, lag.max = 365.24219*24, main='')
# title(main='Stage ACF (8760 hours)', line = 0.5)
# dev.off()

AR_test_forecasts <- list(); AR_models <- list(); params <- list();
# full data - ARIMA(5,1,4)
# full data-100 days - ARIMA(3,1,2)
# 
# Sys.time()
# auto_AR <- auto.arima(dat_train$stageNov, xreg = as.matrix(dat_train[,1:5]), seasonal=T, max.p=36, max.q=36,
#                      max.d = 1, stepwise = T, parallel = F, nmodels=200, test='kpss')  ## -> (3,1,2)
# Sys.time()
ARIMA_orders <- list(c(5,1,1),c(5,1,1),c(5,1,1))#, c(12,1,12)); 
seas <- c(F,T,T)#,F); 
seas_orders <- list(c(0,0,0),c(2,0,1),c(2,1,1))#,c(0,0,0))

dat_stage_ts <- ts(data=dat[(1):N_train,7], frequency=24, class="matrix") 
decomposed <- decompose(dat_stage_ts, type = "additive")
seasonal_component <- decomposed$seasonal
kpss.test(seasonal_component)

dat_stage_ts <- ts(data=dat[(1):N_train,7], frequency=365.24219*24, class="matrix") 
decomposed <- decompose(dat_stage_ts, type = "additive")
seasonal_component <- decomposed$seasonal
kpss.test(seasonal_component)

stl_decomposed <- stl(dat_stage_ts, s.window = "periodic")
seasonal_component <- stl_decomposed$time.series[, "seasonal"]

Sys.time()
AR_models[[1]] <- Arima(dat_train$stageNov, order = c(5,1,1), xreg = as.matrix(dat_train[,1:5]))
Sys.time()
### ARIMA parameters chosen by auto.arima. Seasonal parameters chosen by PACF, spikes at 24 and 48 hours
AR_models[[2]] <- Arima(dat_train$stageNov, order = c(5,1,1), xreg = as.matrix(dat_train[,1:5]),
                        seasonal = list(order=c(2,0,1),period=24))

RMSE_test_nw <- matrix(NA, ncol=length(AR_models), nrow=3); prop_CI_nw <- matrix(NA, ncol=length(AR_models), nrow=3);
colnames(RMSE_test_nw) <- paste0('ARIMAX',params);colnames(prop_CI_nw) <- paste0('ARIMAX',params);
rownames(RMSE_test_nw) <- c('FullDataRMSE', 'Above0RMSE','Below0RMSE');rownames(prop_CI_nw) <- c('FullDataPropPI', 'AboveSeaLevel','BelowSeaLevel')
for(i in 1:length(AR_models)){
  params[[i]] <- paste0(ARIMA_orders[[i]][1],',',ARIMA_orders[[i]][2],',',ARIMA_orders[[i]][3],'-',seas_orders[[i]][1],',',seas_orders[[i]][2],',',seas_orders[[i]][3])
  AR_test_forecasts[[i]] <- forecast(AR_models[[i]], xreg=as.matrix(dat_test[,1:5]))
  
  RMSE_test_nw[1,i] <- sqrt(mean((dat_test$stageNov-AR_test_forecasts[[i]]$mean)^2))
  prop_CI_nw[1,i]   <- length(intersect(which(dat_test$stageNov>AR_test_forecasts[[i]]$lower[,2]),which(dat_test$stageNov<AR_test_forecasts[[i]]$upper[,2])))/N_test

  aboveInd <- which(dat_test$stageNov>=0); aboveN <- length(aboveInd)
  RMSE_test_nw[2,i] <- sqrt(mean((dat_test$stageNov[aboveInd]-AR_test_forecasts[[i]]$mean[aboveInd])^2))
  prop_CI_nw[2,i]   <- length(intersect(which(dat_test$stageNov[aboveInd]>AR_test_forecasts[[i]]$lower[aboveInd,2]),which(dat_test$stageNov[aboveInd]<AR_test_forecasts[[i]]$upper[aboveInd,2])))/aboveN
  
  belowInd <- which(dat_test$stageNov<0); belowN <- length(belowInd)
  RMSE_test_nw[3,i] <- sqrt(mean((dat_test$stageNov[belowInd]-AR_test_forecasts[[i]]$mean[belowInd])^2))
  prop_CI_nw[3,i]   <- length(intersect(which(dat_test$stageNov[belowInd]>AR_test_forecasts[[i]]$lower[belowInd,2]),which(dat_test$stageNov[belowInd]<AR_test_forecasts[[i]]$upper[belowInd,2])))/belowN
}
RMSE_test_nw; prop_CI_nw;

### AR IS MODEL METRICS
AR_model_metrics <- data.frame(cbind(paste0('ARIMAX(',params[[1]],')'), round(accuracy(AR_models[[1]]),6), round(AIC(AR_models[[1]]),2)));
for(i in 2:length(AR_models)){
  metrics <- data.frame(cbind(paste0('ARIMAX(',params[[i]],')'), round(accuracy(AR_models[[i]]),6), round(AIC(AR_models[[i]]),2)));
  AR_model_metrics <- rbind(AR_model_metrics, metrics)
}; colnames(AR_model_metrics) <- c('Model', 'ME', 'RMSE', 'MAE', 'MPE', 'MAPE', 'MASE', 'ACF1', 'AIC')
AR_model_metrics

library(flextable); library(officer); ft <- flextable(as.data.frame(AR_model_metrics)); doc <- read_docx(); doc <- body_add_flextable(doc, value = ft); print(doc, target = "AR_model_metrics.docx")

## DENSE FORECASTING
start_tm <- Sys.time()
print(start_tm); modelNames <-  paste0('ARIMAX',params);
N_models=length(AR_models); N_horizon=96; horizons <- c(1,24,48,72,96)
forecasts <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1))); 
forecast_residuals <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon+1)));
forecasts_intervals <- lapply(1:N_models, function(x) list());
forecast_rmses <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon)));
forecast_maes <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_bounds <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon))); 
forecast_horizon_rmses <- matrix(NA, nrow=N_models, ncol=5); 
rownames(forecast_horizon_rmses) <- modelNames; colnames(forecast_horizon_rmses) <- paste0(horizons,'hr')
max_resids <- rep(NA, N_models);min_resids <- rep(NA, N_models);mean_rmses <- rep(NA, N_models);sd_rmses <- rep(NA, N_models);mean_maes <- rep(NA, N_models);mean_bounds <- rep(NA, N_models);
forecast_long_horizon <- list();forecast_residual_long_horizon <- list();forecasts_interval_long_horizon <- list();
forecast_rmse_long_horizon <- rep(NA,N_models);forecast_mae_long_horizon <- rep(NA,N_models);forecast_bound_long_horizon <- rep(NA,N_models);
inits <- as.ts(dat_train$stageNov[(N_train-23):N_train])
set.seed(96)
for(m in 1:N_models){
  print(paste0('Model: ', params[[m]]))
  AR <- AR_models[[m]]; train_dat <- dat_train
  #restart <- ifelse(m==1, 156, 1)
  # forecasts[[m]] <- cbind(forecasts[[m]],rep(0,96))
  # forecast_residuals[[m]] <- cbind(forecast_residuals[[m]],rep(0,96))
  # last_known <- cbind(dat_test[1:2352,],stage_test[1:2352]); colnames(last_known) <- colnames(train_dat)
  # train_dat <- rbind(train_dat, last_known)
  for(fc in 1:(N_test-N_horizon+1)){
    print(paste0('Forecast Iter: ', m,'-',fc))
    horizon <- fc:(fc+N_horizon-1)
    forecast_m <- forecast(AR, xreg=as.matrix(dat_test[horizon,]), level=0.95)
    forecasts[[m]][,fc] <- forecast_m$mean
    forecasts_intervals[[m]][[fc]] <- cbind(forecast_m$lower, forecast_m$upper)
    resid_m <-  stage_test[horizon] - forecast_m$mean
    forecast_residuals[[m]][,fc] <- resid_m
    forecast_rmses[[m]][fc] <- sqrt(mean(resid_m^2))
    forecast_maes[[m]][fc] <- mean(abs(resid_m))
    forecast_bounds[[m]][fc] <- round(length(intersect(which(stage_test[horizon]>forecast_m$lower),which(stage_test[horizon]<forecast_m$upper)))/N_horizon,4)
    
    last_known <- cbind(dat_test[fc,],stage_test[fc]); colnames(last_known) <- colnames(train_dat)
    train_dat <- rbind(train_dat, last_known)
    AR <- Arima(train_dat$stageNov, order = ARIMA_orders[[m]], seasonal = seas_orders[[m]], xreg = as.matrix(train_dat[,1:5]))
  }
  
  for(h in 1:length(horizons)){
    forecast_horizon_rmses[m,h] <- sqrt(mean(forecast_residuals[[m]][horizons[h],]^2))
  }
  mean_rmses[m] <- round(mean(forecast_rmses[[m]]),6)
  sd_rmses[m] <- round(sd(forecast_rmses[[m]]),6)
  mean_maes[m] <- round(mean(forecast_maes[[m]]),6)
  mean_bounds[m] <- round(mean(forecast_bounds[[m]]),6)
  min_resids[m] <- round(min(forecast_residuals[[m]]),6)
  max_resids[m] <- round(max(forecast_residuals[[m]]),6)

  horizon=1:N_test
  forecast_m <- forecast(AR_models[[m]], xreg=as.matrix(dat_test[horizon,]), level=0.95)
  forecast_long_horizon[[m]] <- forecast_m$mean
  forecasts_interval_long_horizon[[m]] <- cbind(forecast_m$lower, forecast_m$upper)
  resid_m <-  stage_test[horizon] - forecast_m$mean
  forecast_residual_long_horizon[[m]] <- resid_m
  forecast_rmse_long_horizon[m] <- sqrt(mean(resid_m^2))
  forecast_mae_long_horizon[m] <- mean(abs(resid_m))
  forecast_bound_long_horizon[m] <- round(length(intersect(which(stage_test[horizon]>forecast_m$lower),which(stage_test[horizon]<forecast_m$upper)))/N_horizon,4)
}
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

forecast_horizon_rmses
t(t(mean_rmses))
t(t(sd_rmses))
t(t(mean_maes))
t(t(mean_bounds))
t(t(min_resids))
t(t(max_resids))

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



lowTideRMSEs <- c()
highTideBin <- matrix(NA, nrow=96, ncol=2352)
for(t in 1:(2448-96)){
  highTideBin[,t] <- (stage_test[t:(t+95)]>0)
}
for(i in 1:length(forecast_residuals)){
  ht_rmse_sum <- c()
  for(j in 1:2352){ht_rmse_sum <- c(ht_rmse_sum, sqrt(mean((forecast_residuals[[i]][,j][highTideBin[,j]])^2)))}
  print(mean(ht_rmse_sum)*100)
  lt_rmse_sum <- c()
  for(j in 1:2352){
    lt_rmse_sum <- c(lt_rmse_sum, sqrt(mean((forecast_residuals[[i]][,j][!highTideBin[,j]])^2)))
    sqrt(mean((forecast_residuals[[i]][,j][!highTideBin[,j]])^2))
  }
  print(mean(lt_rmse_sum)*100)
}


n_forecast=1; horizon=(n_forecast):(n_forecast+95); i=2#2205 #420
rng=range(c(forecasts[[i]][,n_forecast],stage_test[horizon]))

pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/96-hour Forecasts/',modelNames[i],'-',n_forecast,'.pdf'), width = 9, height = 4)
par(mfrow=c(1,1), mar=c(4,4,2,0.5)); 
plot(stage_test[horizon]~dt_test[horizon], col='darkgrey', pch=19, cex=1.2,
     main=paste0(modelNames[i],' 96-hr Forecast starting ',dt_test[horizon][1]), 
     ylab='Stage (m)',xlab='Datetime',
     ylim=rng)
lines(forecasts[[i]][,n_forecast]~dt_test[horizon], col='blue', lwd=2)
dev.off()

