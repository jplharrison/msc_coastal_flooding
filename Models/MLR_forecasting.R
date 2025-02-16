rm(list=ls())

# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#save.image(file='MLR_scale_forecasting.RData')
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
          end=c(2022,6477.562), frequency=365.24219*24, class="matrix") #end=c(2021,12,28,15)
#dat <- as.data.frame(dat)
tm_yr <- time(dat); start_tm <- min(tm_yr); end_tm <- max(tm_yr); yrs <- end_tm-start_tm; tm_yr <- ((tm_yr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
tm_28 <- time(dat); start_tm <- min(tm_28); end_tm <- max(tm_28); yrs <- (end_tm-start_tm)*365.24219/27.33; tm_28 <- ((tm_28-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#tm_14 <- time(dat); start_tm <- min(tm_14); end_tm <- max(tm_14); yrs <- (end_tm-start_tm)*365.24219/14; tm_14 <- ((tm_14-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#dat <- as.data.frame(cbind(sin(tm_14),cos(tm_14),sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:6] <- c('sin_14','cos_14','sin_28','cos_28','sin_yr','cos_yr')
dat <- as.data.frame(cbind(sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:4] <- c('sin_28','cos_28','sin_yr','cos_yr')
# 32756     9
N <- dim(dat)[1]
N_test <- 2448
N_train <- N-N_test
dat_train <- dat[1:N_train,]
dt_train <- dt[1:N_train]
dat_test <- dat[(N_train+1):N,1:8]
resNov_test <- dat[(N_train+1):N,9]
dt_test <- dt[(N_train+1):N]
var_names <- c('Atmospheric Pressure', 'Local Wind', 'Napa Flow', 'Ocean Wind', 'Residual Water Level')


# OLS Multiple Linear Regression model
linReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind, data = dat_train)
resLinReg <- ts(data=resid(linReg), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")
#stepAIC(linReg)

# Multiple Quadratic Regression
quadReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind+
                I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(oceanWind^2)+
                I(atmPres*gnossWind)+I(atmPres*napaFlow)+I(atmPres*oceanWind)+
                I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+
                I(napaFlow*oceanWind)
              , data = dat_train)
resQuadReg <- ts(data=resid(quadReg), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")
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
cubReg <- lm(formula = resNov ~ atmPres+gnossWind+napaFlow+oceanWind+
     I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(atmPres^3)+
     I(gnossWind^3)+I(napaFlow^3)+I(atmPres*gnossWind)+
     I(atmPres*napaFlow)+I(atmPres*oceanWind)+I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+I(atmPres^2*gnossWind)+
     I(atmPres^2*napaFlow)+I(atmPres^2*oceanWind)+I(gnossWind^2*atmPres)+I(gnossWind^2*napaFlow)+I(gnossWind^2*oceanWind)+
     I(napaFlow^2*atmPres)+I(napaFlow^2*gnossWind)+I(oceanWind^2*atmPres)+I(oceanWind^2*gnossWind)+I(oceanWind^2*napaFlow)+
     I(atmPres*gnossWind*napaFlow)+I(gnossWind*napaFlow*oceanWind), data = dat_train)
resCubReg <- ts(data=resid(cubReg), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")

# OLS Multiple Linear Regression model with Seasonality
linRegSeas <- lm(resNov~sin_28+cos_28+sin_yr+cos_yr+atmPres+gnossWind+napaFlow+oceanWind, data = dat_train)
#stepAIC(linRegSeas)
resLinRegSeas <- ts(data=resid(linRegSeas), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")


# Cubic Regression with Seasonality 
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
                 , data = dat_train)

#stepAIC(cubRegSeas)
cubRegSeas <- lm(resNov ~ sin_28+cos_28+sin_yr+cos_yr+
                   atmPres+gnossWind+napaFlow+oceanWind+
                   I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(oceanWind^2)+
                   I(atmPres^3)+               I(napaFlow^3)+I(oceanWind^3)+
                   I(atmPres*gnossWind)+I(atmPres*napaFlow)+I(atmPres*oceanWind)+
                   I(gnossWind*napaFlow)+                   I(napaFlow*oceanWind)+
                   I(atmPres^2*gnossWind)+I(atmPres^2*napaFlow)+I(atmPres^2*oceanWind)+
                   I(gnossWind^2*atmPres)+                      I(gnossWind^2*oceanWind)+
                   I(napaFlow^2*atmPres)+                       I(napaFlow^2*oceanWind)+
                   I(oceanWind^2*atmPres)+I(oceanWind^2*gnossWind)+
                   I(atmPres*gnossWind*napaFlow)+I(atmPres*gnossWind*oceanWind)+
                                                 I(gnossWind*napaFlow*oceanWind)
                 , data = dat_train)
resCubRegSeas <- ts(data=resid(cubRegSeas), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")

dfs=20
splineAICs <- matrix(ncol=4, nrow=dfs)
for(df in 1:dfs){
  atmPnss <- ns(dat_train$atmPres, df=df)
  localWnss <- ns(dat_train$gnossWind, df=df)
  napaFnss <- ns(dat_train$napaFlow, df=df)
  oceanWnss <- ns(dat_train$oceanWind, df=df)
  splineAICs[df,1] <- AIC(lm(dat_train$resNov~atmPnss))
  splineAICs[df,2] <- AIC(lm(dat_train$resNov~localWnss))
  splineAICs[df,3] <- AIC(lm(dat_train$resNov~napaFnss))
  splineAICs[df,4] <- AIC(lm(dat_train$resNov~oceanWnss))
}
par(mfrow=c(2,2)); for(i in 1:4){
  splineAICs[,i] <- rescale(splineAICs[,i], to=c(0,1))
  plot(splineAICs[,i], type='b', ylab='Relative AIC Proportion', xlab='Number of Degrees of Freedom', main=var_names[i], pch=19); abline(h=0.1, col='lightgrey', lty=2, lwd=1); abline(h=0.2, col='lightgrey', lty=2, lwd=1)}
for(i in 1:4){print(paste0(colnames(dat_train)[i+4], ' DF=',min(which(splineAICs[,i]<0.2))))}
## These are the smallest DFs that capture >90% / >80% of the improvement in AIC of a 20DF spline

m_spline <- lm(resNov~ns(atmPres, df = 2)+ns(gnossWind, df = 3)+ns(napaFlow, df = 6)+ns(oceanWind, df = 2), data=dat_train)
resSpline <- ts(data=resid(m_spline), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")

m_splineSeas <- lm(resNov~sin_28+cos_28+sin_yr+cos_yr+ns(atmPres, df = 2)+ns(gnossWind, df = 3)+ns(napaFlow, df = 6)+ns(oceanWind, df = 2), data=dat_train)
resSplineSeas <- ts(data=resid(m_splineSeas), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")


m_spline2 <- lm(resNov~ns(atmPres, df = 3)+ns(gnossWind, df = 3)+ns(napaFlow, df = 11)+ns(oceanWind, df = 3), data=dat_train)
resSpline2 <- ts(data=resid(m_spline2), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")

m_spline2Seas <- lm(resNov~sin_28+cos_28+sin_yr+cos_yr+ns(atmPres, df = 3)+ns(gnossWind, df = 3)+ns(napaFlow, df = 11)+ns(oceanWind, df = 3), data=dat_train)
resSpline2Seas <- ts(data=resid(m_spline2Seas), start=c(2019, 20), end=c(2022,6477.562-N_test), frequency=365.24219*24, class="matrix")


# Plots of Fitted values for 'Linear', 'Quadratic','Cubic','Linear+Seasonal','Cubic+Seasonal' models
models=list(linReg, quadReg, cubReg, linRegSeas, cubRegSeas, m_spline,m_splineSeas, m_spline2,m_spline2Seas)#,linRegSeas)
modelNames <- c('Linear', 'Quadratic','Cubic','Linear+Seasonal','Cubic+Seasonal', 'Splines(2,3,6,2)', 'Splines(2,3,6,2)+Seas', 'Splines(3,3,11,3)', 'Splines(3,3,11,3)+Seas')
N_models <- length(models); modelCodes <- paste0('M1:',1:N_models)
AICs <- vector(length=length(models));ranks <- vector(length=length(models));BICs <- vector(length=length(models));MSEs <- vector(length=length(models)); RMSEs <- vector(length=length(models)); SSEs <- vector(length=length(models)); minRes <- vector(length=length(models)); maxRes <- vector(length=length(models));
resids=list(resLinReg, resQuadReg, resCubReg, resLinRegSeas, resCubRegSeas, resSpline,resSplineSeas, resSpline2,resSpline2Seas)
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
  MSEs[i] <- round(mean((resids[[i]])^2),4)
  RMSEs[i] <- round(sqrt(mean((resids[[i]])^2)),4)
  minRes[i] <- round(min(resids[[i]]),2)
  maxRes[i] <- round(max(resids[[i]]),2)
  #plot(resids[[i]][1:N_plot]~dt[1:N_plot], type='l', ylim=c(-0.4, 0.3), main=paste0(modelNames[i],'-',round(SSEs[i],2),'-(',minRes[i],',',maxRes[i],')'), xlab='date', ylab='residuals')
  #abline(h=0, lty=2, col='grey')
}
# Residuals QQ plots
par(mfrow=c(2,ceiling(length(models)/2)))
for(i in 1:length(models)){
  qqPlot(resids[[i]], main=modelNames[i], ylab='Residuals', xlab='Normal Quantiles', id=F)
  #qqnorm(resids[[i]], main=modelNames[i])
  #qqline(resids[[i]])
}


modelMetrics <- data.frame(SSEs, MSEs, RMSEs, AICs, BICs, ranks); rownames(modelMetrics) <- modelNames; modelMetrics
library(flextable); library(officer); ft <- flextable(as.data.frame(modelMetrics)); doc <- read_docx(); doc <- body_add_flextable(doc, value = ft); print(doc, target = "MLR_model_metrics.docx")

set.seed(37)
### MLR Forecast
forecasts <- list(); fc_conf <- list();
fc_errors <- list(); forecast_metrics <- list();
N_window <- c(12,24,72); K=length(N_window);
for(i in 1:N_models){
  forecast_metrics[[i]] <- data.frame(Unadjusted=rep(0,4), 
                                      RollingWindow12=rep(0,4), 
                                      PeriodicWindow12=rep(0,4), 
                                      RollingWindow24=rep(0,4), 
                                      PeriodicWindow24=rep(0,4), 
                                      RollingWindow72=rep(0,4), 
                                      PeriodicWindow72=rep(0,4)); 
  rownames(forecast_metrics[[i]]) <- c('ErrorMean', 'ErrorStdDev', 'RMSE', 'PropConfBound')
  
  forecast_m <- predict(models[[i]], newdata = dat_test, interval = "prediction", level=0.95)
  forecasts[[i]] <- forecast_m[,1]
  fc_errors[[i]] <- resNov_test-forecasts[[i]]
  fc_conf[[i]] <- forecast_m[,2:3]
  
  forecast_metrics[[i]][1,1] <- round(mean(fc_errors[[i]]),4)
  forecast_metrics[[i]][2,1] <- round(sd(fc_errors[[i]]),4)
  forecast_metrics[[i]][3,1] <- round(sqrt(mean((fc_errors[[i]])^2)),4)
  forecast_metrics[[i]][4,1] <- round(length(intersect(which(resNov_test>fc_conf[[i]][,1]),which(resNov_test<fc_conf[[i]][,2])))/N_test,4)
}

### MLR + Rolling Window MA Error Correction
forecasts_corrected_errors <- list(); forecasts_corrected <- list(); fc_conf_corrected <- list()
rolling_ma_log <- matrix(NA, nrow=(N_test), ncol=N_models*K); 
for(i in 1:N_models){
  for(k in 1:K){
    resids_window <- resids[[i]][(N_train-N_window[k]+1):N_train]
    resids_ma <- mean(resids_window)
    forecasts_corrected[[(i-1)*K+k]] <- vector(length=N_test)
    fc_conf_corrected[[(i-1)*K+k]] <- matrix(nrow=N_test, ncol=2)
    for(j in 1:N_test){
      rolling_ma_log[j,(i-1)*K+k] <- resids_ma
      fc <- predict(models[[i]], newdata = dat_test[j,], interval = "prediction", level=0.95)
      forecasts_corrected[[(i-1)*K+k]][j] <- fc[,1] + resids_ma
      fc_conf_corrected[[(i-1)*K+k]][j,1:2] <- fc[,2:3]+resids_ma
      resids_window <- c(resids_window[2:N_window[k]], fc)
      resids_ma <- mean(resids_window)
    }
    forecasts_corrected_errors[[(i-1)*K+k]] <- resNov_test-forecasts_corrected[[(i-1)*K+k]]
    forecast_metrics[[i]][1,2*k] <- round(mean(forecasts_corrected_errors[[(i-1)*K+k]]),4)
    forecast_metrics[[i]][2,2*k] <- round(sd(forecasts_corrected_errors[[(i-1)*K+k]]),4)
    forecast_metrics[[i]][3,2*k] <- round(sqrt(mean((forecasts_corrected_errors[[(i-1)*K+k]])^2)),4)
    forecast_metrics[[i]][4,2*k] <- round(length(intersect(which(resNov_test>fc_conf_corrected[[(i-1)*K+k]][,1]),which(resNov_test<fc_conf_corrected[[(i-1)*K+k]][,2])))/N_test,4)
  }  
}

### MLR Forecast + Periodic Window MA Error Correction
forecasts_corrected_errors_fixed <- list(); forecasts_corrected_fixed <- list(); fc_conf_corrected_fixed <- list()
periodic_ma_log <- list(pma12=matrix(NA, nrow=(N_test/N_window[1]), ncol=N_models),
                        pma24=matrix(NA, nrow=(N_test/N_window[2]), ncol=N_models),
                        pma72=matrix(NA, nrow=(N_test/N_window[3]), ncol=N_models))
for(i in 1:N_models){
  for(k in 1:K){
    resids_window <- resids[[i]][(N_train-N_window[k]+1):N_train]
    resids_ma <- mean(resids_window)
    forecasts_corrected_fixed[[(i-1)*K+k]] <- vector(length=N_test)
    fc_conf_corrected_fixed[[(i-1)*K+k]] <- matrix(nrow=N_test, ncol=2)
    for(j in 1:(N_test/N_window[k])){
      periodic_ma_log[[k]][j,i] <- resids_ma
      fc <- predict(models[[i]], newdata = dat_test[((j-1)*N_window[k]+1):(j*N_window[k]),], interval = "prediction", level=0.95)
      forecasts_corrected_fixed[[(i-1)*K+k]][((j-1)*N_window[k]+1):(j*N_window[k])] <- fc[,1] + resids_ma
      fc_conf_corrected_fixed[[(i-1)*K+k]][((j-1)*N_window[k]+1):(j*N_window[k]),1:2] <- fc[,2:3]+resids_ma
      resids_window <- fc[,1]
      resids_ma <- mean(resids_window)
    }
    forecasts_corrected_errors_fixed[[(i-1)*K+k]] <- resNov_test-forecasts_corrected_fixed[[(i-1)*K+k]]
    forecast_metrics[[i]][1,(1+2*k)] <- round(mean(forecasts_corrected_errors_fixed[[(i-1)*K+k]]),4)
    forecast_metrics[[i]][2,(1+2*k)] <- round(sd(forecasts_corrected_errors_fixed[[(i-1)*K+k]]),4)
    forecast_metrics[[i]][3,(1+2*k)] <- round(sqrt(mean((forecasts_corrected_errors_fixed[[(i-1)*K+k]])^2)),4)
    forecast_metrics[[i]][4,(1+2*k)] <- round(length(intersect(which(resNov_test>fc_conf_corrected_fixed[[(i-1)*K+k]][,1]),which(resNov_test<fc_conf_corrected_fixed[[(i-1)*K+k]][,2])))/N_test,4)
  }
}
### Add pred bounds
forecastMetrics <- list(RMSE= as.data.frame(matrix(NA, ncol=N_models+1, nrow=1+K*2)), PropBound= as.data.frame(matrix(NA, ncol=N_models+1, nrow=1+K*2))); fcMs <- c('RMSE','PropBound')
for(j in 1:2){
  colnames(forecastMetrics[[j]]) <- c(fcMs[j],modelNames); forecastMetrics[[j]][,1] <- c('Unadjusted','RollingWindow12','PeriodicWindow12','RollingWindow24','PeriodicWindow24','RollingWindow72','PeriodicWindow72')
  for(i in 1:N_models){
    forecastMetrics[[j]][,(i+1)] <- unlist(as.vector(forecast_metrics[[i]][(2+j),]))
    print(modelNames[i])
    print(forecast_metrics[[i]])
  }
}
forecastMetrics

### Using OOS test data, 
##### the linear model performs the best, both in RMSE and PropBound
##### the non-seasonal models outperform the models with the seasonal component
##### the splines models are outperformed by the polynomial models
##### the error correction improves the non-seasonal polynomial models
##### but error correction worsens the performance of the seasonal polynomial models as well as all of the spline models

ft <- flextable(as.data.frame(forecastMetrics)); doc <- read_docx(); doc <- body_add_flextable(doc, value = ft); print(doc, target = "MLR_forecast_metrics_.docx")

par(mfrow=c(3,3)); k=3; i=1;
for(i in 1:3){
  for(k in 1:K){
    plot(rep(periodic_ma_log[[k]][,i], each=N_window[k])~dt_test, type='l', lwd=2, ylab='Error Correction', main=paste0(modelNames[i], ' - ', N_window[k], ' hour window'))
    lines(rolling_ma_log[,(i-1)*K+k]~dt_test, type='l', lwd=2, lty=2, col='red')
  }
}
par(mfrow=c(3,3))
## MLR 100 Day Forecast

for(i in 1:N_models){
  for(k in 1:K){
    plot(dat$resNov[(N-2500):(N)]~dt[(N-2500):(N)], type='l', lwd=4, col='darkgrey', 
         main=paste0('Test Forecast: ',modelNames[i]), ylab='resNov', xlab='date')
    lines(forecasts[[(i-1)*K+k]]~dt_test, col=1, lty=1, lwd=2)
  }
}


### MLR ~3 week forecast
plot(dat$resNov[(N-2500):(N-1800)]~dt[(N-2500):(N-1800)], type='l', lwd=3, col='darkgrey', 
     main='Test Forecast', ylab='resNov', xlab='date')
abline(v=dt_test[1], col='lightgrey', lty=2,lwd=2)
for(i in 1:N_models){lines(forecasts[[i]][1:600]~dt_test[1:600], col=i, lty=2, lwd=2)}
legend('topleft', legend=modelNames, lwd=2, lty=2, col=1:N_models, bty='n')

## MLR Forecast Errors
par(mfrow=c(1,2))
for(i in 1:N_models){
  plot(fc_errors[[i]]~dt_test, type='l', lwd=2, ylim=c(-0.25,0.25),#ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0(modelNames[i], ' Forecast Error'), ylab='Forecast Error', xlab='Date')
  legend('topleft', legend=c(paste0('Error Mean: ', round(mean(fc_errors[[i]]),3)), paste0('Error Std Dev: ', round(sd(fc_errors[[i]]),3))), bty='n')
  plot(forecasts_corrected_errors[[i]]~dt_test, type='l', lwd=2, ylim=c(-0.25,0.25),#ylim=c(min(unlist(forecasts_corrected_errors)), max(unlist(forecasts_corrected_errors))),
       main=paste0(modelNames[i], ' Corrected Forecast Error'), ylab='Forecast Error', xlab='Date')
  legend('topleft', legend=c(paste0('Error Mean: ', round(mean(forecasts_corrected_errors[[i]]),3)), paste0('Error Std Dev: ', round(sd(forecasts_corrected_errors[[i]]),3))), bty='n')
}

## Confidence interval for PolyReg is really tight - bounds capture very few of the observed points. Look at Conf/Var adjustment
## 4 day Forecast: Unadjusted vs Rolling Window Error correction vs Periodic Error Correction
forecast_period <- 7*24; legend_cex=0.9;
modelInd <- c(1,6,8)
par(mfrow=c(length(modelInd),3))
for(i in modelInd){#N_models){
  ylims <- c(min(c(unlist(fc_conf), unlist(fc_conf_corrected), unlist(fc_conf_corrected_fixed))), max(c(unlist(fc_conf), unlist(fc_conf_corrected), unlist(fc_conf_corrected_fixed)))) #ylim=c(min(c(unlist(forecasts),resNov_test)), max(c(unlist(forecasts),resNov_test)))
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=4, 
       ylim=ylims, xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=,c(-0.25,0.25)
       main=paste0(modelNames[i], ' Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test[1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=4, col='darkgrey')
  lines(forecasts[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=3, col='blue')
  lines(fc_conf[[i]][1:forecast_period,1]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  lines(fc_conf[[i]][1:forecast_period,2]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  abline(h=0, lty=2, col='lightgrey', lwd=2)
  legend('topleft', legend=c(paste0('RMSE: ', forecast_metrics[[i]][3,1]),
                             paste0('Prop Conf. Bound: ', forecast_metrics[[i]][4,1])), bty='n', pch=19, cex=legend_cex)
  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=4, 
       ylim=ylims, xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=,c(-0.25,0.25)
       main=paste0(modelNames[i], ' Rolling Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test[1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=4, col='darkgrey')
  lines(forecasts_corrected[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=3, col='blue')
  lines(fc_conf_corrected[[i]][1:forecast_period,1]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  lines(fc_conf_corrected[[i]][1:forecast_period,2]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  abline(h=0, lty=2, col='lightgrey', lwd=2)
  legend('topleft', legend=c(paste0('RMSE: ', forecast_metrics[[i]][3,2]),
                             paste0('Prop Conf. Bound: ', forecast_metrics[[i]][4,2])), bty='n', pch=19, cex=legend_cex)

  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=4, 
       ylim=ylims, xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=,c(-0.25,0.25)
       main=paste0(modelNames[i], ' Periodic Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test[1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=4, col='darkgrey')
  lines(forecasts_corrected_fixed[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=3, col='blue')
  lines(fc_conf_corrected_fixed[[i]][1:forecast_period,1]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  lines(fc_conf_corrected_fixed[[i]][1:forecast_period,2]~dt_test[1:forecast_period], type='l', lwd=2, col='red')
  abline(h=0, lty=2, col='lightgrey', lwd=2)
  legend('topleft', legend=c(paste0('RMSE: ', forecast_metrics[[i]][3,3]),
                             paste0('Prop Conf. Bound: ', forecast_metrics[[i]][4,3])), bty='n', pch=19, cex=legend_cex)
}

## 4 day Forecast Error: Unadjusted vs Rolling Window vs Periodic Window Error correction
modelInd <- c(1,5,6)
par(mfrow=c(length(modelInd),3))
for(i in modelInd){#N_models){
  plot(resids[[i]][(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.12,0.12), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0(modelNames[i], ' Forecast Error'), ylab='Forecast Error', xlab='Date')
  lines(fc_errors[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2)
  abline(h=0, lty=2, col='lightgrey', lwd=2)
   
  plot(resids[[i]][(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.12,0.12), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Rolling Correction Forecast Error'), ylab='Forecast Error', xlab='Date')
  lines(forecasts_corrected_errors[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2)
  abline(h=0, lty=2, col='lightgrey', lwd=2)
  #legend('topleft', legend=c(paste0('Error Mean: ', round(mean(forecasts_corrected_errors[[i]][1:forecast_period]),3)), 
  #                           paste0('Error Std Dev: ', round(sd(forecasts_corrected_errors[[i]][1:forecast_period]),3)), 
  #                           paste0('RMSE: ', round(sqrt(sum((forecasts_corrected_errors[[i]][1:forecast_period])^2)/forecast_period),3))), bty='n')
  
  
  plot(resids[[i]][(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.12,0.12), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Periodic Correction Forecast Error'), ylab='Forecast Error', xlab='Date')
  lines(forecasts_corrected_errors_fixed[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2)
  abline(h=0, lty=2, col='lightgrey', lwd=2)
  #legend('topleft', legend=c(paste0('Error Mean: ', round(mean(forecasts_corrected_errors_fixed[[i]][1:forecast_period]),3)), 
  #                           paste0('Error Std Dev: ', round(sd(forecasts_corrected_errors_fixed[[i]][1:forecast_period]),3)), 
  #                           paste0('RMSE: ', round(sqrt(sum((forecasts_corrected_errors_fixed[[i]][1:forecast_period])^2)/forecast_period),3))), bty='n')
}





## 4 day Forecasts: Unadjusted vs Rolling Window vs Periodic Window Error correction
par(mfrow=c(N_models,3)); forecast_period <- 4*24
for(i in 1:N_models){
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.2), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0(modelNames[i], ' Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.2), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Rolling Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts_corrected[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.2), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Periodic Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts_corrected_fixed[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
}


## 4 day Forecasts for Linear+Seasonal: Unadjusted vs Rolling Window vs Periodic Window Error correction
par(mfrow=c(2,3)); forecast_period <- 100*24
for(i in c(1,4)){
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.45), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0(modelNames[i], ' Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.45), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Rolling Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts_corrected[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
  
  plot(dat_train$resNov[(N_train-N_window):N_train]~dt_train[(N_train-N_window):N_train], type='l', lwd=2, ylim=c(-0.2,0.45), xlim=c(dt_train[N_train-12], dt_test[forecast_period]), col='brown', #ylim=c(min(unlist(fc_errors)), max(unlist(fc_errors))),
       main=paste0('Periodic Correction Forecast'), ylab='resNov', xlab='Date')
  lines(resNov_test~dt_test, lwd=2)
  lines(forecasts_corrected_fixed[[i]][1:forecast_period]~dt_test[1:forecast_period], type='l', lwd=2, lty=2, col='blue')
}





