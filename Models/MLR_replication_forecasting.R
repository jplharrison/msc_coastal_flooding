rm(list=ls())

# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='MLR_forecasting.RData')
#load(file='MLR_forecasting.RData')
#load(file='MLR_forecasting.RData')

#Libraries
library(car)
library(forecast)
library(ggplot2)
library(gridExtra)
library(MASS)
library(scales)
library(splines)
library(tseries)
library(tidyverse)
library(urca)
library(vars)
library(zoo)
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

# Create time series object
dt <- dat_nov[,1]
dat <- ts(data=dat_nov[2:6], start=c(2019, 20),
          end=c(2022,6477.562), frequency=365.24219*24, class="matrix") #end=c(2021,12,28,15)
dat_stage <- ts(data=dat_nov[7:8], start=c(2019, 20),
          end=c(2022,6477.562), frequency=365.24219*24, class="matrix")
#dat <- as.data.frame(dat)
# tm_yr <- time(dat); start_tm <- min(tm_yr); end_tm <- max(tm_yr); yrs <- end_tm-start_tm; tm_yr <- ((tm_yr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
# tm_28 <- time(dat); start_tm <- min(tm_28); end_tm <- max(tm_28); yrs <- (end_tm-start_tm)*365.24219/27.33; tm_28 <- ((tm_28-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#tm_14 <- time(dat); start_tm <- min(tm_14); end_tm <- max(tm_14); yrs <- (end_tm-start_tm)*365.24219/14; tm_14 <- ((tm_14-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#dat <- as.data.frame(cbind(sin(tm_14),cos(tm_14),sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:6] <- c('sin_14','cos_14','sin_28','cos_28','sin_yr','cos_yr')
#dat <- as.data.frame(cbind(sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:4] <- c('sin_28','cos_28','sin_yr','cos_yr')
dat <- as.data.frame(dat)
# 32756     9
N <- dim(dat)[1]
N_test <- 2448
N_train <- N-N_test
dat_train <- dat[1:N_train,]
dt_train <- dt[1:N_train]
dat_test <- dat[(N_train+1):N,1:4]
dat_stage_train <- dat_stage[1:N_train,]
dat_stage_test <- dat_stage[(N_train+1):N,]
resNov_test <- dat[(N_train+1):N,5]
dt_test <- dt[(N_train+1):N]
var_names <- c('Atmospheric Pressure', 'Local Wind', 'Napa Flow', 'Ocean Wind', 'Residual Water Level')


# OLS Multiple Linear Regression model
linReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind, data = dat_train)
resLinReg <- ts(data=resid(linReg), start=c(2019, 20),
                end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")
#stepAIC(linReg)

# Multiple Quadratic Regression
quadReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind+
                I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(oceanWind^2)+
                I(atmPres*gnossWind)+I(atmPres*napaFlow)+I(atmPres*oceanWind)+
                I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+
                I(napaFlow*oceanWind)
              , data = dat_train)
resQuadReg <- ts(data=resid(quadReg), start=c(2019, 20),
                 end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")
#stepAIC(quadReg)

# Multiple Quadratic Regression
cubReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind+
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
             , data = dat_train)
#stepAIC(cubReg)
# cubReg <- lm(formula = resNov ~ atmPres+gnossWind+napaFlow+oceanWind+
#                I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(atmPres^3)+
#                I(gnossWind^3)+I(napaFlow^3)+I(atmPres*gnossWind)+
#                I(atmPres*napaFlow)+I(atmPres*oceanWind)+I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+I(atmPres^2*gnossWind)+
#                I(atmPres^2*napaFlow)+I(atmPres^2*oceanWind)+I(gnossWind^2*atmPres)+I(gnossWind^2*napaFlow)+I(gnossWind^2*oceanWind)+
#                I(napaFlow^2*atmPres)+I(napaFlow^2*gnossWind)+I(oceanWind^2*atmPres)+I(oceanWind^2*gnossWind)+I(oceanWind^2*napaFlow)+
#                I(atmPres*gnossWind*napaFlow)+I(gnossWind*napaFlow*oceanWind), data = dat_train)
resCubReg <- ts(data=resid(cubReg), start=c(2019, 20),
                end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")


# Plots of Fitted values for 'Linear', 'Quadratic','Cubic','Linear+Seasonal','Cubic+Seasonal' models
models=list(linReg, quadReg, cubReg)#,linRegSeas)
modelNames <- c('Linear', 'Quadratic','Cubic')
N_models <- length(models); modelCodes <- paste0('M1:',1:N_models)
AICs <- vector(length=length(models));ranks <- vector(length=length(models));BICs <- vector(length=length(models));MSEs <- vector(length=length(models)); RMSEs <- vector(length=length(models)); SSEs <- vector(length=length(models)); minRes <- vector(length=length(models)); maxRes <- vector(length=length(models));
resids=list(resLinReg, resQuadReg, resCubReg)
preds <- list(); N_plot <- N_train
par(mfrow=c(2,length(models)))
for(i in 1:length(models)){
  preds[[i]] <- predict(models[[i]])
  #plot(dat$resNov[1:N_plot]~dt[1:N_plot], type='l', lwd=4, col='darkgrey', ylim=c(-.2,0.65), main=modelNames[i], ylab='resNov', xlab='date')
  #points(preds[[i]][1:N_plot]~dt[1:N_plot], type='l', col='blue',lty=1, lwd=2)
  AICs[i] <- round(AIC(models[[i]]),2)
  BICs[i] <- round(BIC(models[[i]]),2)
}
# Plots of Residuals for 5 models
for(i in 1:length(resids)){
  ranks[i] <- models[[i]]$rank
  SSEs[i] <- round(sum((resids[[i]])^2),2)
  MSEs[i] <- round(mean((resids[[i]])^2),6)
  RMSEs[i] <- round(sqrt(mean((resids[[i]])^2)),6)
  minRes[i] <- round(min(resids[[i]]),2)
  maxRes[i] <- round(max(resids[[i]]),2)
  #plot(resids[[i]][1:N_plot]~dt[1:N_plot], type='l', ylim=c(-0.4, 0.3), main=paste0(modelNames[i],'-',round(SSEs[i],2),'-(',minRes[i],',',maxRes[i],')'), xlab='date', ylab='residuals')
  #abline(h=0, lty=2, col='grey')
}
# Residuals QQ plots
par(mfrow=c(1,ceiling(length(models)/1)), mar=c(4,4,2,1))
for(i in 1:length(models)){
  qqPlot(resids[[i]], main=paste0(modelNames[i], ' Model'), ylab='Residuals', xlab='Normal Quantiles', id=F)
  #qqnorm(resids[[i]], main=modelNames[i])
  #qqline(resids[[i]])
}


modelMetrics <- data.frame(SSEs, MSEs, RMSEs, AICs, BICs, ranks); rownames(modelMetrics) <- modelNames; modelMetrics
library(flextable); library(officer); ft <- flextable(as.data.frame(modelMetrics)); doc <- read_docx(); doc <- body_add_flextable(doc, value = ft); print(doc, target = "MLR_model_metrics.docx")

set.seed(37)
### MLR Forecast
start_tm <- Sys.time()
print(start_tm)
N_horizon=96; horizons <- c(1,12,24,48,72,96); N_correction=12;
forecasts <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon-N_correction))); 
forecast_residuals <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon-N_correction)));
forecasts_intervals <- lapply(1:N_models, function(x) list());
forecast_rmses <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction)));
forecast_maes <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction))); 
forecast_bounds <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction))); 
forecast_horizon_rmses <- matrix(NA, nrow=N_models, ncol=length(horizons)); rownames(forecast_horizon_rmses) <- modelNames; colnames(forecast_horizon_rmses) <- paste0(horizons,'hr')
max_resids <- rep(NA, N_models);min_resids <- rep(NA, N_models);mean_rmses <- rep(NA, N_models);sd_rmses <- rep(NA, N_models);mean_maes <- rep(NA, N_models);mean_bounds <- rep(NA, N_models);
forecast_long_horizon <- list();forecast_residual_long_horizon <- list();forecasts_interval_long_horizon <- list();
forecast_rmse_long_horizon <- rep(NA,N_models);forecast_mae_long_horizon <- rep(NA,N_models);forecast_bound_long_horizon <- rep(NA,N_models);

forecast_adj_residuals <- lapply(1:N_models, function(x) matrix(NA, nrow = N_horizon, ncol = (N_test - N_horizon-N_correction)));
forecast_adj_rmses <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction)));
forecast_adj_maes <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction))); 
forecast_adj_bounds <- lapply(1:N_models, function(x) rep(NA, (N_test - N_horizon-N_correction))); 
max_resids_adj <- rep(NA, N_models);min_resids_adj <- rep(NA, N_models);mean_rmses_adj <- rep(NA, N_models);sd_rmses_adj <- rep(NA, N_models);mean_maes_adj <- rep(NA, N_models);mean_bounds_adj <- rep(NA, N_models);
forecast_adj_horizon_rmses <- matrix(NA, nrow=N_models, ncol=length(horizons)); rownames(forecast_adj_horizon_rmses) <- modelNames; colnames(forecast_adj_horizon_rmses) <- paste0(horizons,'hr')

set.seed(96)
for(m in 1:N_models){
  
  for(fc in (N_correction+1):(N_test-N_horizon)){
    correction_horizon <- (fc-N_correction):(fc-1)
    horizon <- fc:(fc+N_horizon-1)
    correction_m <- predict(models[[m]], newdata=dat_test[correction_horizon,], interval = "prediction", level=0.95)[,1]
    correction_adj <- mean(resNov_test[correction_horizon]-correction_m)
    forecast_m <- predict(models[[m]], newdata=dat_test[horizon,], interval = "prediction", level=0.95)
    forecast_stage <- forecast_m[,1] + dat_stage_test[horizon, 2]
    forecast_stage_adj <- forecast_stage+correction_adj
    forecasts[[m]][,(fc-N_correction)] <- forecast_stage
    forecast_stage_interval <- forecast_m[,2:3] + dat_stage_test[horizon, 2]
    forecast_stage_interval_adj <- forecast_m[,2:3] + correction_adj + dat_stage_test[horizon, 2]
    forecasts_intervals[[m]][[fc]] <- forecast_stage_interval
    resid_m <-  dat_stage_test[horizon,1] - forecast_stage
    resid_m_adj <-  dat_stage_test[horizon,1] - forecast_stage_adj
    forecast_residuals[[m]][,(fc-N_correction)] <- resid_m
    forecast_rmses[[m]][(fc-N_correction)] <- sqrt(mean(resid_m^2))
    forecast_maes[[m]][(fc-N_correction)] <- mean(abs(resid_m))
    forecast_bounds[[m]][(fc-N_correction)] <- round(length(intersect(which(dat_stage_test[horizon,1]>forecast_stage_interval[,1]),which(dat_stage_test[horizon,1]<forecast_stage_interval[,2])))/N_horizon,6)
    
    forecast_adj_residuals[[m]][,(fc-N_correction)] <- resid_m_adj
    forecast_adj_rmses[[m]][(fc-N_correction)] <- sqrt(mean(resid_m_adj^2))
    forecast_adj_maes[[m]][(fc-N_correction)] <- mean(abs(resid_m_adj))
    forecast_adj_bounds[[m]][(fc-N_correction)] <- round(length(intersect(which(dat_stage_test[horizon,1]>forecast_stage_interval_adj[,1]),which(dat_stage_test[horizon,1]<forecast_stage_interval_adj[,2])))/N_horizon,6)
  }
  
  
  for(h in 1:length(horizons)){
    forecast_horizon_rmses[m,h] <- sqrt(mean(forecast_residuals[[m]][horizons[h],]^2))
    forecast_adj_horizon_rmses[m,h] <- sqrt(mean(forecast_adj_residuals[[m]][horizons[h],]^2))
  }
  mean_rmses[m] <- round(mean(forecast_rmses[[m]]),6)
  sd_rmses[m] <- round(sd(forecast_rmses[[m]]),6)
  mean_maes[m] <- round(mean(forecast_maes[[m]]),6)
  mean_bounds[m] <- round(mean(forecast_bounds[[m]]),6)
  min_resids[m] <- round(min(forecast_residuals[[m]]),6)
  max_resids[m] <- round(max(forecast_residuals[[m]]),6)
  
  mean_rmses_adj[m] <- round(mean(forecast_adj_rmses[[m]]),6)
  sd_rmses_adj[m] <- round(sd(forecast_adj_rmses[[m]]),6)
  mean_maes_adj[m] <- round(mean(forecast_adj_maes[[m]]),6)
  mean_bounds_adj[m] <- round(mean(forecast_adj_bounds[[m]]),6)
  min_resids_adj[m] <- round(min(forecast_adj_residuals[[m]]),6)
  max_resids_adj[m] <- round(max(forecast_adj_residuals[[m]]),6)
  
  horizon=1:N_test
  forecast_m <- predict(models[[m]], newdata=dat_test[horizon,], interval = "prediction", level=0.95)
  forecast_long_stage <- forecast_m[,1] + dat_stage_test[horizon, 2]
  forecast_long_horizon[[m]] <- forecast_long_stage
  forecast_long_stage_interval <- forecast_m[,2:3] + dat_stage_test[horizon, 2]
  forecasts_interval_long_horizon[[m]] <- forecast_long_stage_interval
  resid_m <-  dat_stage_test[horizon,1] - forecast_long_stage
  forecast_residual_long_horizon[[m]] <- resid_m
  forecast_rmse_long_horizon[m] <- sqrt(mean(resid_m^2))
  forecast_mae_long_horizon[m] <- mean(abs(resid_m))
  forecast_bound_long_horizon[m] <- round(length(intersect(which(dat_stage_test[horizon,1]>forecast_long_stage_interval[,1]),which(dat_stage_test[horizon,1]<forecast_long_stage_interval[,2])))/N_horizon,6)
  
}
end_tm <- Sys.time()
print(end_tm)
print(paste0('Duration: ', end_tm-start_tm))

forecast_horizon_rmses
forecast_adj_horizon_rmses
mean_rmses*100
mean_rmses_adj*100

t(t(round(100*mean_rmses,1)))
t(t(round(100*sd_rmses,1)))
t(t(round(100*mean_maes,1)))
t(t(round(100*mean_bounds,1)))
t(t(round(100*min_resids,1)))
t(t(round(100*max_resids,1)))


t(t(mean_rmses_adj))
t(t(sd_rmses_adj))
t(t(mean_maes_adj))
t(t(mean_bounds_adj))
t(t(min_resids_adj))
t(t(max_resids_adj))
