rm(list=ls())


# Data Directory
#setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='ARMIAX_forecasts102.RData') # save.image(file='ARMAX_train.RData')
#load(file='ARMIAX_forecasts102.RData')
#load(file = 'ARMAX_train.RData')

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
dat_nov <- dat_nov[,1:6]
dat_pet <- merge(X, y_pet, by='datetime')
dat_pet <- dat_pet[,1:6]
dat_row <- merge(X, y_row, by='datetime')
colnames(dat_nov) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resNov')
colnames(dat_pet) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resPet')
colnames(dat_row) <- c('datetime', 'atmPres', 'gnossWind','napaFlow','oceanWind', 'resRow')
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
dt <- dat_nov[1:N_train,1]
dt_test <- dat_nov[(1+N_train):N,1]
dat <- ts(data=dat_nov[1:N_train,2:6], start=c(2019, 20),
          end=c(2022,(6477.562-N_test)), frequency=365.24219*24, class="matrix") 
dat_test <- ts(data=dat_nov[(1+N_train):N,2:6], start=c(2022,(6477.562-N_test+1)),
          end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix") 
dat <- as.data.frame(dat); dat_test <- as.data.frame(dat_test);

adf.test(dat$resNov, alternative = "stationary")
#adf.test(dat$resNov, alternative = "explosive")
pa <- pacf(dat$resNov, lag.max = 48)
a <- acf(dat$resNov, lag.max = 24000)

AR_test_forecasts <- list(); AR_models <- list(); params <- list();
# full data - ARIMA(5,1,4)
# full data-100 days - ARIMA(3,1,2)
# 
# auto_AR <- auto.arima(dat$resNov, xreg = as.matrix(dat[,1:4]), seasonal=T, max.p=48, max.q=48,
#                       max.order=80, max.d = 12, stepwise = F, parallel = F)  ## -> (3,1,2)
Sys.time()
AR_models[[1]] <- Arima(dat$resNov, order = c(3,1,2), xreg = as.matrix(dat[,1:4]))
Sys.time()
AR_models[[2]] <- Arima(dat$resNov, order = c(31,1,2), xreg = as.matrix(dat[,1:4]))
Sys.time()
AR_models[[3]] <- Arima(dat$resNov, order = c(31,1,4), xreg = as.matrix(dat[,1:4]))
Sys.time()
AR_models[[4]] <- Arima(dat$resNov, order = c(31,1,12), xreg = as.matrix(dat[,1:4]))
Sys.time()
AR_models[[5]] <- Arima(dat$resNov, order = c(31,1,24), xreg = as.matrix(dat[,1:4]))
Sys.time()
for(i in 1:length(AR_models)){
  params[[i]] <- paste0(length(AR_models[[i]]$model$phi),',',length(AR_models[[i]]$model$Delta),',',length(which(abs(AR_models[[i]]$model$theta)>0)))
  AR_test_forecasts[[i]] <- forecast(AR_models[[i]], xreg=as.matrix(dat_test[,1:4]))
}

### AR IS MODEL METRICS
AR_model_metrics <- data.frame(cbind(paste0('ARIMAX(',params[[1]],')'), round(accuracy(AR_models[[1]]),12), round(AIC(AR_models[[1]]),2)));
for(i in 2:length(AR_models)){
  metrics <- data.frame(cbind(paste0('ARIMAX(',params[[i]],')'), round(accuracy(AR_models[[i]]),12), round(AIC(AR_models[[i]]),2)));
  AR_model_metrics <- rbind(AR_model_metrics, metrics)
}; colnames(AR_model_metrics) <- c('Model', 'ME', 'RMSE', 'MAE', 'MPE', 'MAPE', 'MASE', 'ACF1', 'AIC')
AR_model_metrics

library(flextable); library(officer); ft <- flextable(as.data.frame(AR_model_metrics)); doc <- read_docx(); doc <- body_add_flextable(doc, value = ft); print(doc, target = "AR_model_metrics.docx")


### AR FORECASTING
forecast_periods <- c(12,24,72); K=length(forecast_periods);
rolling_forecast <- list(); rolling_forecast_lo_95 <- list(); rolling_forecast_hi_95 <- list(); RMSE_test=matrix(NA, ncol=K, nrow=length(AR_models)); prop_CI=matrix(NA, ncol=K, nrow=length(AR_models))
train_dat <- dat
#Arima(dat$resNov, order = c(3,1,2), xreg = as.matrix(dat[,1:4]))
for(i in 1:length(AR_models)){
  print(paste0('Model: ', params[[i]]))
  for(k in 1:K){
    print(paste0('Forecast period: ', forecast_periods[k]))
    rolling_forecast[[(i-1)*K+k]] <- vector(); rolling_forecast_lo_95[[(i-1)*K+k]] <- vector(); rolling_forecast_hi_95[[(i-1)*K+k]] <- vector();
    AR <- AR_models[[i]]; train_dat <- dat
    for(j in 1:(N_test/forecast_periods[k])){
      print(j)
      print(dim(train_dat))
      fc <- forecast(AR, xreg=as.matrix(dat_test[(((j-1)*forecast_periods[k])+1):(j*forecast_periods[k]),1:4]), h = forecast_periods[k])
      rolling_forecast[[(i-1)*K+k]] <- c(rolling_forecast[[(i-1)*K+k]], fc$mean)
      rolling_forecast_lo_95[[(i-1)*K+k]] <- c(rolling_forecast_lo_95[[(i-1)*K+k]], fc$lower[,2])
      rolling_forecast_hi_95[[(i-1)*K+k]] <- c(rolling_forecast_hi_95[[(i-1)*K+k]], fc$upper[,2])
      train_dat <- rbind(train_dat, dat_test[(((j-1)*forecast_periods[k])+1):(j*forecast_periods[k]),])
      AR <- Arima(train_dat$resNov, order = as.numeric(unlist(strsplit(params[[i]], ','))), xreg = as.matrix(train_dat[,1:4]))
    }
    RMSE_test[i,k] <- sqrt(mean((dat_test$resNov-rolling_forecast[[(i-1)*K+k]])^2))
    prop_CI[i,k]   <- length(intersect(which(dat_test$resNov>rolling_forecast_lo_95[[(i-1)*K+k]]),which(dat_test$resNov<rolling_forecast_hi_95[[(i-1)*K+k]])))/N_test
  }
}  
# Proportion of test obs inside 95% confidence interval
length(intersect(which(dat_test$resNov>AR_test_forecasts[[2]]$lower[,2]),which(dat_test$resNov<AR_test_forecasts[[2]]$upper[,2])))/N_test
length(intersect(which(dat_test$resNov>rolling_forecast_lo_95[[1]]),which(dat_test$resNov<rolling_forecast_hi_95[[1]])))/N_test
length(intersect(which(dat_test$resNov>rolling_forecast_lo_95[[2]]),which(dat_test$resNov<rolling_forecast_hi_95[[2]])))/N_test
length(intersect(which(dat_test$resNov>rolling_forecast_lo_95[[3]]),which(dat_test$resNov<rolling_forecast_hi_95[[3]])))/N_test


# Rolling Forecasts AR(3,1,2)
portrait=F; par(mfrow=c(ifelse(portrait,3,1),ifelse(portrait,1,3))); plot_n <- ifelse(portrait,24*102,24*9); ablines=T; abline_period=ifelse(portrait, 24*7, 24); i=1; 
for(k in 1:3){
  plot(dat_test$resNov[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3, col='darkgrey', 
       main=paste0('ARIMAX(',params[[i]],') - ',forecast_periods[k], ' Hour Prediction Window'), xlab='Datetime', ylab='resNov',
       xlim=c(dt[(N_train-5)],dt_test[plot_n]), ylim=c(-.5,.5))
  lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
  lines(rolling_forecast[[k]][1:plot_n]~dt_test[1:plot_n], col='blue', lty=1,lwd=2)
  lines(rolling_forecast_lo_95[[k]][1:plot_n]~dt_test[1:plot_n], col='red', lty=1,lwd=1)
  lines(rolling_forecast_hi_95[[k]][1:plot_n]~dt_test[1:plot_n], col='red', lty=1,lwd=1)
  text(x=dt[(N_train-(plot_n/27.2))], y=ifelse(portrait,-0.3,-0.36), cex=ifelse(portrait,0.8,1), adj=c(0,0), labels = c(paste0('Test RMSE: ',round(RMSE_test[k],4))))
  text(x=dt[(N_train-(plot_n/27.2))], y=-0.4, cex=ifelse(portrait,0.8,1), adj=c(0,0), labels = c(paste0('Prop. in Pred. Bound: ', round(prop_CI[k],4))))
  legend('topleft', legend=c('Observed','Forecast','Prediction Bound'), lty=1, bg='white', x.intersp=ifelse(portrait,0.4,0.8), y.intersp=ifelse(portrait,0.9,1),
         col=c('darkgrey','blue','red'), lwd=c(3,2,1), cex=ifelse(portrait,0.8,1), bty='n', inset=0)
  if(ablines){abline(v=dt_test[seq(from=9, by=abline_period, to=plot_n)], lty=2, col='grey')
              abline(v=dt_test[seq(from=9, by=abline_period, to=(abline_period*ifelse(portrait,1,4)))], lty=2, col='white')
    }
}


# ARIMAX312_predictions_residuals
par(mfrow=c(2,2))
# 100 hours testing AR
plot(dat_test$resNov[1:100]~dt_test[1:100], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params[[1]],') Prediction - 100 hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[100]), ylim=c(-.15,.15))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(AR_test_forecasts[[1]]$mean[1:100]~dt_test[1:100], col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Train Obs Tail', 'Test Obs', 'Forecast'), lty=c(1,1,2), col=c('brown', 'darkgrey','blue'), lwd=2)
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours testing AR_models[[2]]
plot(dat_test$resNov[1:100]~dt_test[1:100], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params[[2]],') Prediction - 100 hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[100]), ylim=c(-.15,.15))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(AR_test_forecasts[[2]]$mean[1:100]~dt_test[1:100], col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Train Obs Tail', 'Test Obs', 'Forecast'), lty=c(1,1,2), col=c('brown', 'darkgrey','blue'), lwd=2)
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours testing AR3
plot(dat_test$resNov[1:100]~dt_test[1:100], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params[[3]],') Prediction - 100 hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[100]), ylim=c(-.15,.15))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(AR_test_forecasts[[3]]$mean[1:100]~dt_test[1:100], col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Train Obs Tail', 'Test Obs', 'Forecast'), lty=c(1,1,2), col=c('brown', 'darkgrey','blue'), lwd=2)
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours testing AR4
plot(dat_test$resNov[1:100]~dt_test[1:100], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params[[4]],') Prediction - 100 hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[100]), ylim=c(-.15,.15))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(AR_test_forecasts[[4]]$mean[1:100]~dt_test[1:100], col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Train Obs Tail', 'Test Obs', 'Forecast'), lty=c(1,1,2), col=c('brown', 'darkgrey','blue'), lwd=2)
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')

par(mfrow=c(2,2))
# 100 hours resid AR
plot((dat_test$resNov[1:100]-AR_test_forecasts[[1]]$mean[1:100])~dt_test[1:100], 
     type='l', lwd=3, xlab='Datetime',main=paste0('ARIMAX(',params[[1]],') Test Residual - 100 hours'), ylab='Test Residual', 
     ylim=c(-.05, 0.25), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(0:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours resid AR_models[[2]]
plot((dat_test$resNov[1:100]-AR_test_forecasts[[2]]$mean[1:100])~dt_test[1:100], 
     type='l', lwd=3, xlab='Datetime',main=paste0('ARIMAX(',params[[2]],') Test Residual - 100 hours'), ylab='Test Residual', 
     ylim=c(-.05, 0.25), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(0:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours resid AR3
plot((dat_test$resNov[1:100]-AR_test_forecasts[[3]]$mean[1:100])~dt_test[1:100], 
     type='l', lwd=3, xlab='Datetime',main=paste0('ARIMAX(',params[[3]],') Test Residual - 100 hours'), ylab='Test Residual', 
     ylim=c(-.05, 0.25), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(0:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
# 100 hours resid AR4
plot((dat_test$resNov[1:100]-AR_test_forecasts[[4]]$mean[1:100])~dt_test[1:100], 
     type='l', lwd=3, xlab='Datetime',main=paste0('ARIMAX(',params[[4]],') Test Residual - 100 hours'), ylab='Test Residual', 
     ylim=c(-.05, 0.25), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(0:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')

par(mfrow=c(1,1))
# 100 hours resid AR1-4
plot((dat_test$resNov[1:100]-AR_test_forecasts[[1]]$mean[1:100])~dt_test[1:100], lty=1,
     type='l', lwd=3, xlab='Datetime',main=paste0('ARIMAX(',params[[1]],') Test Residual - 100 hours'), ylab='Test Residual', 
     ylim=c(-.3, 0.3), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(0:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
for(i in 2:4){lines((dat_test$resNov[1:100]-AR_test_forecasts[[i]]$mean[1:100])~dt_test[1:100], lwd=3, lty=2, col=i)}
legend('topleft', legend=c(params[[1]], params[[2]], params[[3]], params[[4]]), lty=2, col=1:4, lwd=3, bty='n')

rolling_forecast2 <- vector(); rolling_forecast_lo_952 <- vector(); rolling_forecast_hi_952 <- vector();p
train_dat <- dat
#Arima(dat$resNov, order = c(3,1,2), xreg = as.matrix(dat[,1:4]))
for(j in 1:100){
  print(j)
  print(dim(train_dat))
  fc <- forecast(AR2, xreg=as.matrix(dat_test[(((j-1)*24)+1):(j*24),1:4]), h = 24)
  rolling_forecast2 <- c(rolling_forecast, fc$mean)
  rolling_forecast_lo_952 <- c(rolling_forecast_lo_952, fc$lower[,2])
  rolling_forecast_hi_952 <- c(rolling_forecast_hi_952, fc$upper[,2])
  train_dat <- rbind(train_dat, dat_test[(((j-1)*24)+1):(j*24),])
  AR2 <- Arima(train_dat$resNov, order = c(31,1,12), xreg = as.matrix(train_dat[,1:4]))
}


i=2; 
### 
plot_n <- 2400; ablines=T; abline_period=28*24; plot_static_prediction=F;
par(mfrow=c(1,2))
# 284 hours testing AR(3,1,2)
if(plot_static_prediction){par(mfrow=c(1,3))
plot(dat_test$resNov[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params[[2]],') Prediction - ',plot_n,' hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[plot_n]), ylim=c(-.5,.5))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(AR_test_forecasts[[2]]$mean[1:plot_n]~dt_test[1:plot_n], col='blue', lty=2,lwd=2)
lines(AR_test_forecasts[[2]]$lower[1:plot_n,2]~dt_test[1:plot_n], col='red', lty=2,lwd=1)
lines(AR_test_forecasts[[2]]$upper[1:plot_n,2]~dt_test[1:plot_n], col='red', lty=2,lwd=1)
if(ablines){abline(v=dt_test[seq(from=1, by=abline_period, to=plot_n)], lty=2, col='grey')}}






### 284 day analysis: rapid drop in Ocean Wind and Rapid increase in Local Wind
par(mfrow=c(1,5))
for(k in 1){
  plot(dat_test$resNov[1:plot_n]~dt_test[1:plot_n], type='l', lwd=4, col='darkgrey', 
       main=paste0('ARIMAX(',params[[2]],') Rolling ',forecast_periods[k], ' Hourly Prediction - ',plot_n,' hours'), xlab='Datetime', ylab='resNov',
       xlim=c(dt[(N_train-5)],dt_test[plot_n]), ylim=c(-.5,.5))
  lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=4, col='brown')
  lines(rolling_forecast[[k]][1:plot_n]~dt_test[1:plot_n], col='blue', lty=1,lwd=2.5)
  lines(rolling_forecast_lo_95[[k]][1:plot_n]~dt_test[1:plot_n], col='red', lty=1,lwd=1)
  lines(rolling_forecast_hi_95[[k]][1:plot_n]~dt_test[1:plot_n], col='red', lty=1,lwd=1)
  abline(v=dt_test[seq(from=1, by=24, to=plot_n)], lty=2, col='grey')
}
plot(dat_test$atmPres[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3); abline(v=dt_test[seq(from=1, by=24, to=plot_n)], lty=2, col='grey')
plot(dat_test$gnossWind[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3); abline(v=dt_test[seq(from=1, by=24, to=plot_n)], lty=2, col='grey')
plot(dat_test$napaFlow[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3); abline(v=dt_test[seq(from=1, by=24, to=plot_n)], lty=2, col='grey')
plot(dat_test$oceanWind[1:plot_n]~dt_test[1:plot_n], type='l', lwd=3); abline(v=dt_test[seq(from=1, by=24, to=plot_n)], lty=2, col='grey')




acf2(resCubRegSeas, max.lag = 8670)

# OLS Multiple Linear Regression model
linReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind, data = dat)
resLinReg <- ts(data=resid(linReg) , start=c(2019, 20), end=c(2022,6477.562), frequency=365.24219*24, class="matrix")

cubRegSeas <- lm(resNov ~ sin_28+cos_28+sin_yr+cos_yr+
                   atmPres+gnossWind+napaFlow+oceanWind+
                   I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(oceanWind^2)+
                   I(atmPres^3)+I(gnossWind^3)+I(napaFlow^3)+I(oceanWind^3)+
                   I(atmPres*gnossWind)+I(atmPres*napaFlow)+I(atmPres*oceanWind)+
                   I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+I(napaFlow*oceanWind)+
                   I(atmPres^2*gnossWind)+I(atmPres^2*napaFlow)+I(atmPres^2*oceanWind)+
                   I(gnossWind^2*atmPres)+I(gnossWind^2*napaFlow)+I(gnossWind^2*oceanWind)+
                   I(napaFlow^2*atmPres)+I(napaFlow^2*gnossWind)+I(napaFlow^2*oceanWind)+
                   I(oceanWind^2*atmPres)+I(oceanWind^2*gnossWind)+I(oceanWind^2*napaFlow)+
                   I(atmPres*gnossWind*napaFlow)+I(atmPres*gnossWind*oceanWind)+
                   I(atmPres*napaFlow*oceanWind)+I(gnossWind*napaFlow*oceanWind)
                 , data = dat)
resCubRegSeas <- ts(data=resid(cubRegSeas) , start=c(2019, 20), end=c(2022,6477.562), frequency=365.24219*24, class="matrix")

acf2(resLinReg, max.lag = 36, ylim=c(-.1,.1))
acf2(resCubRegSeas, max.lag = 8670)

ord=2;

AR2Pred <- predict(linReg2, newxreg = cbind(dat$atmPres, dat$gnossWind, dat$napaFlow, dat$oceanWind))

par(mfrow=c(1,1))
plot(dat$resNov~dt, pch=19, col='darkgrey', main=paste0('ARX',ord,' Prediction'))
lines(AR2Pred$pred~dt, type='l', ylim=c(-.2,.7), col='blue')

# full data - ARIMA(5,1,4)
# full data-100 days - ARIMA(3,1,2)
AR_auto <- auto.arima(dat$resNov, xreg = as.matrix(dat[,1:4]))
params <- paste0(length(AR_auto$model$phi),',',length(AR_auto$model$Delta),',',length(AR_auto$model$theta))
ARforecast <- forecast(AR_auto, xreg = as.matrix(dat[,1:4]))
ARforecast_test <- forecast(AR_auto, xreg=as.matrix(dat_test[,1:4]))
accuracy(AR_auto)

par(mfrow=c(1,1))
# Training + Testing
plot(dat$resNov~dt, pch=19, col='darkgrey', main=paste0('ARIMAX(',params,') Prediction'), xlim=c(min(dt),max(dt_test)))
lines(ARforecast$mean~dt, col='red')
lines(ARforecast_test$mean~dt_test, col='blue')

# ARIMAX312_predictions_residuals
par(mfrow=c(2,2))
# 100 hours testing
plot(dat_test$resNov[1:100]~dt_test[1:100], type='l', lwd=3, col='darkgrey', 
     main=paste0('ARIMAX(',params,') Prediction - 100 hours'), xlab='Datetime', ylab='resNov',
     xlim=c(dt[(N_train-5)],dt_test[100]), ylim=c(-.2,.3))
lines(dat$resNov[(N_train-5):N_train]~dt[(N_train-5):N_train], lty=1,lwd=2, col='brown')
lines(ARforecast_test$mean[1:100]~dt_test[1:100], col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Train Obs Tail', 'Test Obs', 'Forecast'), lty=c(1,1,2), col=c('brown', 'darkgrey','blue'), lwd=2)
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
#abline(v=dt_test[1:100], lty=2, col='grey')
# Testing
plot(dat_test$resNov~dt_test, type='l', lwd=3, col='darkgrey', xlab='Datetime',ylab='resNov',main=paste0('ARIMAX(',params,') Prediction - 100 days (2400 hours)'), xlim=c(min(dt_test),max(dt_test)), ylim=c(-.2,.3))
lines(ARforecast_test$mean~dt_test, col='blue', lty=2,lwd=2)
legend('topleft', legend=c('Test Obs', 'Forecast'), lty=c(1,2), col=c('darkgrey','blue'),lwd=2)
abline(v=dt_test[seq(from=1, by=7*24, to=2400)], lty=2, col='grey')
# 100 hours resid
plot((dat_test$resNov[1:100]-ARforecast_test$mean[1:100])~dt_test[1:100], 
     type='l', lwd=3, xlab='Datetime',main='Test Residual - 100 hours', ylab='Test Residual', 
     ylim=c(0, 0.45), xlim=c(dt[(N_train-5)],dt_test[100]))
abline(h=0.1*(1:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=24, to=2400)], lty=2, col='grey')
#legend('topleft', legend=c('Test Residual', '24 hour interval'), lty=c(1,2), col=c('black', 'grey'),lwd=2)
# 100 days resid
plot((dat_test$resNov[1:2400]-ARforecast_test$mean[1:2400])~dt_test[1:2400], type='l', lwd=3, xlab='Datetime',main='Test Residual - 100 days (2400 hours)', ylab='Test Residual', ylim=c(0, 0.45))
abline(h=0.1*(1:4), lty=2, col='grey')
abline(v=dt_test[seq(from=1, by=7*24, to=2400)], lty=2, col='grey')
#legend('topleft', legend=c('Test Residual', '1 week interval'), lty=c(1,2), col=c('black', 'grey'),lwd=2)



