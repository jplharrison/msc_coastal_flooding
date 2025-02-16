rm(list=ls())


# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='gam1.RData')
#load(file = 'gam1.RData')
#save.image(file='gamDenseForecast.RData')
#load(file='gamDenseForecast.RData')
#Libraries
library(forecast)
library(ggplot2)
library(gratia)
library(gridExtra)
library(mgcv)
library(nlme)
library(scales)
library(tidyverse)
library(tseries)

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

N <- dim(X)[1]; N_lags=24; 
dat_lags <- matrix(NA, nrow=N, ncol=N_lags)
for(i in 1:N_lags){dat_lags[(1+i):(N),i]=dat_nov$resNov[(1):(N-i)]}
colnames(dat_lags) <- paste0('lag_',1:N_lags)

dat <- na.omit(cbind(dat_nov%>%select(datetime, atmPres, gnossWind, napaFlow, oceanWind, predNov, stageNov), dat_lags))
dt <- dat$datetime; stage <- dat$stageNov; dat <- dat[,-c(1,7)]


N <- dim(dat)[1]; N_test <- 2448; N_train <- N-N_test

dt_train <- dt[1:N_train]
dt_test <- dt[(1+N_train):N]

dat_train <- as.data.frame(ts(data=dat[1:N_train,], start=c(2019, 20+N_lags),
                              end=c(2022,(4008.562)), frequency=365.24219*24, class="matrix"))
dat_test <- as.data.frame(ts(data=dat[(1+N_train):N,], start=c(2022,(4008.562+1)),
                             end=c(2022,(4008.562+2448)), frequency=365.24219*24, class="matrix"))
stage_train <- stage[1:N_train]
stage_test <- stage[(1+N_train):N]

S <- diag(1)*1000
penParams <- list(m1=NULL,
                  m2=list(lag_1=list(S),lag_2=list(S),lag_3=list(S)),
                  m3=list(lag_1=list(S),lag_2=list(S),lag_3=list(S),
                          lag_4=list(S),lag_5=list(S),lag_6=list(S),
                          lag_7=list(S),lag_8=list(S),lag_9=list(S),
                          lag_10=list(S),lag_11=list(S),lag_12=list(S),
                          lag_13=list(S),lag_14=list(S),lag_15=list(S),
                          lag_16=list(S),lag_17=list(S),lag_18=list(S),
                          lag_19=list(S),lag_20=list(S),lag_21=list(S),
                          lag_22=list(S),lag_23=list(S),lag_24=list(S)))


lag_vector <- c(0,3,24)
gam_models <- list(); gam_summaries <- list(); gam_ISresids <- list(); gam_forecasts <- list(); gam_OOSresids <- list(); gam_OOSinterval <- list();
for(l in 1:length(lag_vector)){
  print(paste0('Lag: ', lag_vector[l], ' - ', Sys.time()))
  # Create time series object
  N <- dim(X)[1]; N_lags=lag_vector[l]; 
  if(N_lags>0){
    train_dat <- dat_train %>% select(atmPres, gnossWind, napaFlow, oceanWind, predNov, paste0('lag_',1:N_lags))
  } else{
    train_dat <- dat_train%>%select(atmPres, gnossWind, napaFlow, oceanWind, predNov)
  }
  
  N <- dim(dat)[1]; N_test <- 2448; N_train <- dim(train_dat)
  
  if(N_lags>0){
  #lag_terms <- paste(paste0("s(", colnames(dat_lags)[1:N_lags], ", bs='ps')"), collapse = " + ")
    lag_terms <- paste(paste0(colnames(dat_lags)[1:N_lags]), collapse = " + ")
    formula_string <- paste("stage_train ~ s(atmPres) + s(gnossWind) + s(napaFlow) + s(oceanWind) + s(predNov) +", lag_terms)
  }else{formula_string <- "stage_train ~ s(atmPres) + s(gnossWind) + s(napaFlow) + s(oceanWind) + s(predNov)"}
  S=diag(1)
  gam_models[[l]] <- gam(as.formula(formula_string),
                        data = train_dat,
                        family = scat(), paraPen = penParams[[l]])  # s(x,k=10)
  gam_summaries[[l]] <- summary(gam_models[[l]] )
  
  # IS Fitted Values
  gam_ISresids[[l]] <- stageNov_train-gam_models[[l]]$fitted.values
  
  ### OOS Prediction and Residuals
  pred <- predict(gam_models[[l]], newdata = dat_test, se.fit=T)
  sigma <- sqrt(sum(residuals(gam_models[[l]])^2) / df.residual(gam_models[[l]]))
  gam_forecasts[[l]] <- pred$fit
  gam_OOSinterval[[l]] <- cbind(pred$fit-1.96*sqrt(pred$se.fit^2+sigma^2), pred$fit+1.96*sqrt(pred$se.fit^2+sigma^2))
  gam_OOSresids[[l]] <- stageNov_test-pred$fit
  
  print(paste0('End - ', Sys.time()))
}  



## DENSE FORECAST
modelNames <-  c('GAM0', 'GAM3', 'GAM24')
N_models=length(gam_models); N_horizon=96; horizons <- c(1,24,48,72,96)
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
N_lags <- c(0,3,24);
start_tm <- Sys.time(); print(start_tm); set.seed(96)
## For each GAM model
for(m in 1){
  N <- dim(X)[1]
  
  dat <- dat_nov%>%select(atmPres, gnossWind, napaFlow, oceanWind, predNov)
  dt <- dat_nov$datetime; stage <- dat_nov$stageNov
  
  N <- dim(dat)[1]; N_test <- 2448; N_train <- N-N_test
  dt_test <- dt[(1+N_train):N]
  test_dat <- as.data.frame(ts(data=dat[(1+N_train):N,], start=c(2022,(6477.562-N_test+1)),
                               end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix"))
  stage_test <- stage[(1+N_train):N]
  
  ## Iterate through 2353 columns
  for(fc in 1:(N_test-N_horizon+1)){
    print(paste0(m,'-',fc))
    ## set horizon window of 96 hours
    horizon <- fc:(fc+N_horizon-1)
    ## iterate through 96 hours of forecast
    X_fc <- test_dat[horizon,]
    ## forecast horizon
    forecast_m <- predict(gam_models[[m]], newdata = X_fc, se.fit=T)
    forecasts[[m]][,fc] <- forecast_m$fit
    ## calculate residuals and metrics for each 96 hour forecast
    resid_m <-  stage_test[horizon] - forecast_m$fit
    forecast_residuals[[m]][,fc] <- resid_m
    
    t_value <- qt(0.975,df=N_train-1)
    sigm <- sd(resid_m)
    forecast_lower <- forecast_m$fit-(t_value+0.2) *sigm
    forecast_upper <- forecast_m$fit+(t_value+0.2) *sigm
    forecast_bounds[[m]] <- cbind(forecast_lower, forecast_upper)
    
    forecast_rmses[[m]][fc] <- sqrt(mean(resid_m^2))
    forecast_maes[[m]][fc] <- mean(abs(resid_m))
    forecast_bounds[[m]][fc] <- round(length(intersect(which(stage_test[horizon]>forecast_lower),which(stage_test[horizon]<forecast_upper)))/N_horizon,4)
  }
  ## calculate aggregate metrics for model
  for(h in 1:length(horizons)){
    forecast_horizon_rmses[m,h] <- sqrt(mean(forecast_residuals[[m]][horizons[h],]^2))
  }
  mean_rmses[m] <- round(mean(forecast_rmses[[m]]),6)
  sd_rmses[m] <- round(sd(forecast_rmses[[m]]),6)
  mean_maes[m] <- round(mean(forecast_maes[[m]]),6)
  mean_bounds <- round(mean(forecast_bounds[[m]]),6)
  min_resids[m] <- round(min(forecast_residuals[[m]]),6)
  max_resids[m] <- round(max(forecast_residuals[[m]]),6)
}
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

start_tm <- Sys.time(); print(start_tm); set.seed(96)
## For each GAM model
for(m in 2:N_models){
  N <- dim(dat_nov)[1]
  dat_lags <- matrix(NA, nrow=N, ncol=N_lags[m])
  for(i in 1:N_lags[m]){
    dat_lags[(1+i):(N),i]=dat_nov$stageNov[(1):(N-i)]
  }
  colnames(dat_lags) <- paste0('lag_',1:N_lags[m]); #dat_lags <- rbind(matrix(NA, nrow=N_lags[m], ncol=N_lags[m]), dat_lags)
  
  dat <- na.omit(cbind(dat_nov%>%select(datetime, atmPres, gnossWind, napaFlow, oceanWind, predNov, stageNov), dat_lags))
  dt <- dat$datetime; stage <- dat$stageNov; dat <- dat[,-c(1,7)]
  
  N <- dim(dat)[1]; N_test <- 2448; N_train <- N-N_test
  dt_test <- dt[(1+N_train):N]
  test_dat <- as.data.frame(ts(data=dat[(1+N_train):N,], start=c(2022,(6477.562-N_test+1)),
                               end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix"))
  stage_test <- stage[(1+N_train):N]
  
  ## Iterate through 2353 columns
  for(fc in 1:(N_test-N_horizon+1)){
    ## set horizon window of 96 hours
    horizon <- fc:(fc+N_horizon-1)
    ## initialise lag responses at start of forecast window
    lag_window <- test_dat[fc,(1+5):(N_lags[m]+5)];
    ## iterate through 96 hours of forecast
    for(i in horizon){
      print(paste0(m,'-',fc,'-',i))
      ## update lagged values for hour of interest + generate dmatrix of one row
      test_dat[i,(1+5):(N_lags[m]+5)] <- lag_window
      X_i <- test_dat[i,]
      ## forecast hour of interest
      forecast_m <- predict(gam_models[[m]], newdata = X_i, se.fit=T)
      forecasts[[m]][(i-fc+1),fc] <- forecast_m$fit
      ## update lagged values with most recent forecast
      lag_window=unlist(c(forecast_m$fit, lag_window))[1:(N_lags[m])];
    }
    ## calculate residuals and metrics for each 96 hour forecast
    resid_m <-  stage_test[horizon] - forecasts[[m]][,fc]
    forecast_residuals[[m]][,fc] <- resid_m
    
    t_value <- qt(0.975,df=N_train-1)
    sigm <- sd(resid_m)
    forecast_m <- forecasts[[m]][,fc]
    forecast_lower <- forecast_m-(t_value+0.2) *sigm
    forecast_upper <- forecast_m+(t_value+0.2) *sigm
    forecast_bounds[[m]] <- cbind(forecast_lower, forecast_upper)
    
    forecast_rmses[[m]][fc] <- sqrt(mean(resid_m^2))
    forecast_maes[[m]][fc] <- mean(abs(resid_m))
    forecast_bounds[[m]][fc] <- round(length(intersect(which(stage_test[horizon]>forecast_lower),which(stage_test[horizon]<forecast_upper)))/N_horizon,4)
    
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

for(m in 2:3){
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


t(t(colMeans(forecast_hightide_rmses)))
t(t(colMeans(forecast_hightide_maes)))
t(forecast_hightide_rmses[horizons,])

sqrt(mean((svm_models[[2]]$residuals)^2))

forecast_horizon_rmses
cbind(mean_rmses, sd_rmses, mean_maes)
t(t(mean_rmses))
t(t(sd_rmses))
t(t(mean_maes))
t(t(mean_bounds))
t(t(min_resids))
t(t(max_resids))












OOSinterval[[l]] <- cbind(gam_forecasts[[l]]$fit-1.96*sqrt(gam_forecasts[[l]]$se.fit^2+sigma^2), gam_forecasts[[l]]$fit+1.96*sqrt(gam_forecasts[[l]]$se.fit^2+sigma^2))
# }

dt_test <- dt[(1+N_train):N]
stageNov_test <- stage[(1+N_train):N]

modelMetrics <- matrix(NA, nrow=4, ncol=length(lag_vector)); 
rownames(modelMetrics) <- c('RMSE','PropPredBound','MinResid','MaxResid'); 
colnames(modelMetrics) <- paste0('GAM',1:length(gam_models),'-',  lag_vector)
for(l in 1:length(lag_vector)){
  modelMetrics[1,l] <- sqrt(mean(gam_OOSresids[[l]]^2))
  modelMetrics[2,l] <- sum((((stageNov_test>=gam_OOSinterval[[l]][,1]) & (stageNov_test<=gam_OOSinterval[[l]][,2]))))/length(stageNov_test)
  modelMetrics[3,l] <- min(gam_OOSresids[[l]])
  modelMetrics[4,l] <- max(gam_OOSresids[[l]])
}; modelMetrics

par(mfrow=c(length(lag_vector),5))
for(l in 1:length(lag_vector)){
  #gam.check(gam_models[[l]])
}


par(mfrow=c(2, ceiling(length(lag_vector)/2)))
for(l in 1:length(lag_vector)){
  plot(stageNov_test~dt_test, type='l', lwd=4, col='darkgrey', ylim=range(c(unlist(gam_forecasts), stageNov_test)),
     main=paste0('Forecast: GAM',l, ', ',lag_vector[l],' hr lag'), ylab='stageNov', xlab='date')
  lines(gam_forecasts[[l]]$fit~dt_test, col=1, lty=1, lwd=2)
  #lines(gam_OOSinterval[[l]][,1]~dt_test, col='red', lty=1, lwd=2)
  #lines(gam_OOSinterval[[l]][,2]~dt_test, col='red', lty=1, lwd=2)
}


par(mfrow=c(2, ceiling(length(lag_vector)/2))); plot_n=72
for(l in 1:length(lag_vector)){
  plot(stageNov_test[1:plot_n]~dt_test[1:plot_n], type='l', lwd=4, col='darkgrey', ylim=range(c(unlist(gam_forecasts), stageNov_test)),
       main=paste0('Forecast: GAM',l, ', ',lag_vector[l],' hr lag'), ylab='stageNov', xlab='date')
  lines(gam_forecasts[[l]]$fit[1:plot_n]~dt_test[1:plot_n], col=1, lty=1, lwd=2)
  lines(gam_OOSinterval[[l]][1:plot_n,1]~dt_test[1:plot_n], col='red', lty=1, lwd=2)
  lines(gam_OOSinterval[[l]][1:plot_n,2]~dt_test[1:plot_n], col='red', lty=1, lwd=2)
}

par(mfrow=c(2, ceiling(length(lag_vector)/2)));
for(l in 1:length(lag_vector)){
  plot(gam_OOSresids[[l]]~dt_test, col=2, lty=1, lwd=2, type='l', ylim=range(c(unlist(gam_OOSresids))),
       main=paste0('OOS Residuals: GAM',l, ', ',lag_vector[l],' hr lag'))
}



post <- gam.mh(gam_models[[1]], ns=10000, burn=2000)


par(mfrow=c(4,3))
for(i in 2:10){hist(post$bs[,i], main=colnames(post$bs)[i])}
for(i in 38:46){hist(post$bs[,i], main=colnames(post$bs)[i])}


draw(gam_models[[1]])

for(l in 1:length(gam_models)){appraise(gam_models[[l]])}

pdf('GAM_appraise0.pdf', height=8, width=8)
appraise(gam_models[[1]])
dev.off()
pdf('GAM_appraise3.pdf', height=8, width=8)
appraise(gam_models[[2]])
dev.off()
pdf('GAM_appraise24.pdf', height=8, width=8)
appraise(gam_models[[3]])
dev.off()
pdf('GAM_appraise4.pdf', height=8, width=8)
appraise(gam_models[[4]])
dev.off()
pdf('GAM_appraise5.pdf', height=8, width=8)
appraise(gam_models[[5]])
dev.off()



