rm(list=ls())

# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")

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
#for(i in 2:5){X[,i] <- rescale(X[,i], to=c(0,1))}
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
          end=c(2022,6477), frequency=365.25*24, class="matrix") #end=c(2021,12,28,15)
#dat <- as.data.frame(dat)
tm_yr <- time(dat); start_tm <- min(tm_yr); end_tm <- max(tm_yr); yrs <- end_tm-start_tm; tm_yr <- ((tm_yr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
tm_28 <- time(dat); start_tm <- min(tm_28); end_tm <- max(tm_28); yrs <- (end_tm-start_tm)*365.25/27.33; tm_28 <- ((tm_28-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#tm_14 <- time(dat); start_tm <- min(tm_14); end_tm <- max(tm_14); yrs <- (end_tm-start_tm)*365.25/14; tm_14 <- ((tm_14-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#dat <- as.data.frame(cbind(sin(tm_14),cos(tm_14),sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:6] <- c('sin_14','cos_14','sin_28','cos_28','sin_yr','cos_yr')
dat <- as.data.frame(cbind(sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:4] <- c('sin_28','cos_28','sin_yr','cos_yr')

par(mfrow=c(3,3))
for(i in 1:9){
  plot(dat[,i]~dt, type='l', main=colnames(dat)[i], ylab=colnames(dat)[i])
}

# OLS Multiple Linear Regression model
linReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind, data = dat)
resLinReg <- ts(data=resid(linReg) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

# Multiple Quadratic Regression
quadReg <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind+
                I(atmPres^2)+I(gnossWind^2)+I(napaFlow^2)+I(oceanWind^2)+
                I(atmPres*gnossWind)+I(atmPres*napaFlow)+I(atmPres*oceanWind)+
                I(gnossWind*napaFlow)+I(gnossWind*oceanWind)+
                I(napaFlow*oceanWind)
                , data = dat)
resQuadReg <- ts(data=resid(quadReg) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

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
                 , data = dat)
resCubReg <- ts(data=resid(cubReg) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

# OLS Multiple Linear Regression model with Seasonality
linRegSeas <- lm(resNov~sin_28+cos_28+sin_yr+cos_yr+atmPres+gnossWind+napaFlow+oceanWind, data = dat)
resLinRegSeas <- ts(data=resid(linRegSeas) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")


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
                   I(atmPres*napaFlow*oceanWind)++I(gnossWind*napaFlow*oceanWind)
                 , data = dat)
resCubRegSeas <- ts(data=resid(cubRegSeas) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

dfs=20
splineAICs <- matrix(ncol=4, nrow=dfs)
for(df in 1:dfs){
  atmPnss <- ns(dat$atmPres, df=df)
  localWnss <- ns(dat$gnossWind, df=df)
  napaFnss <- ns(dat$napaFlow, df=df)
  oceanWnss <- ns(dat$oceanWind, df=df)
  splineAICs[df,1] <- AIC(lm(dat$resNov~atmPnss))
  splineAICs[df,2] <- AIC(lm(dat$resNov~localWnss))
  splineAICs[df,3] <- AIC(lm(dat$resNov~napaFnss))
  splineAICs[df,4] <- AIC(lm(dat$resNov~oceanWnss))
}
par(mfrow=c(2,2)); for(i in 1:4){
  splineAICs[,i] <- rescale(splineAICs[,i], to=c(0,1))
  plot(splineAICs[,i], type='b', ylab='AIC', main=colnames(dat)[i+4]); abline(h=0.2, col='lightgrey', lty=2)}
for(i in 1:4){print(paste0(colnames(dat)[i+4], ' DF=',min(which(splineAICs[,i]<0.2))))}
## These are the smallest DFs that capture >80% of the improvement in AIC of a 20DF spline

m_spline <- lm(resNov~ns(atmPres, df = 2)+ns(gnossWind, df = 3)+ns(napaFlow, df = 6)+ns(oceanWind, df = 2), data=dat)
resSpline <- ts(data=resid(m_spline) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

m_splineSeas <- lm(resNov~sin_28+cos_28+sin_yr+cos_yr+ns(atmPres, df = 2)+ns(gnossWind, df = 3)+ns(napaFlow, df = 6)+ns(oceanWind, df = 2), data=dat)
resSplineSeas <- ts(data=resid(m_spline) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

# Plots of Fitted values for 'Linear', 'Quadratic','Cubic','Linear+Seasonal','Cubic+Seasonal' models
models=list(linReg, quadReg, cubReg, linRegSeas, cubRegSeas, m_spline, m_splineSeas)#,linRegSeas)
modelNames <- c('Linear', 'Quadratic','Cubic','Linear+Seasonal','Cubic+Seasonal', 'Splines', 'Splines+Seasonal')#,'Linear+seasonal')
AICs <- vector(length=length(models));BICs <- vector(length=length(models));MSEs <- vector(length=length(models)); RMSEs <- vector(length=length(models)); SSEs <- vector(length=length(models)); minRes <- vector(length=length(models)); maxRes <- vector(length=length(models));
resids=list(resLinReg, resQuadReg, resCubReg, resLinRegSeas, resCubRegSeas, resSpline, resSplineSeas)#, resLinRegSeas)
preds <- list();
par(mfrow=c(2,length(models)))
for(i in 1:length(models)){
  preds[[i]] <- predict(models[[i]])
  plot(dat$resNov[1:2000]~dt[1:2000], pch=19, col='darkgrey', ylim=c(-.2,0.65), main=modelNames[i], ylab='resNov', xlab='date')
  points(preds[[i]][1:2000]~dt[1:2000], type='l', col='blue',lty=1, lwd=2)
  AICs[i] <- AIC(models[[i]])
  BICs[i] <- BIC(models[[i]])
}
# Plots of Residuals for 5 models
for(i in 1:length(resids)){
  SSEs[i] <- round(sum((resids[[i]])^2),4)
  MSEs[i] <- round(mean((resids[[i]])^2),4)
  RMSEs[i] <- round(sqrt(mean((resids[[i]])^2)),4)
  minRes[i] <- round(min(resids[[i]]),2)
  maxRes[i] <- round(max(resids[[i]]),2)
  plot(resids[[i]][1:2000]~dt[1:2000], type='l', ylim=c(-0.4, 0.2), main=paste0(modelNames[i],'-',round(SSEs[i],2),'-(',minRes[i],',',maxRes[i],')'), xlab='date', ylab='residuals')
  abline(h=0, lty=2, col='grey')
}
# Residuals QQ plots
par(mfrow=c(2,ceiling(length(models)/2)))
for(i in 1:length(models)){
  qqPlot(resids[[i]], main=modelNames[i], ylab='Residuals', xlab='Normal Quantiles', id=F)
  #qqnorm(resids[[i]], main=modelNames[i])
  #qqline(resids[[i]])
}

modelMetrics <- rbind(SSEs, MSEs, RMSEs, AICs, BICs); colnames(modelMetrics) <- modelNames; modelMetrics

dfs=1
par(mfrow=c((dfs+1),1))
ns1 <- ns(dat$atmPres, df=dfs)
plot(dat$atmPres~dt, type='l', lwd=2)
for(i in 1:dfs){plot(ns1[,i]~dt, lwd=2, type='l')}
m_ns1 <- lm(dat$atmPres~ns1); fc_ns1 <- predict(m_ns1, newdata = ns1)
par(mfrow=c(1,1)); 
plot(dat$atmPres~dt, type='l', lwd=3);  lines(fc_ns1~dt, col='red')
plot(dat$atmPres-min(dat$atmPres)~dt, type='l', lwd=3, ylim=c(-10,45))
for(i in 1:dfs){lines(m_ns1$coefficients[i+1]*ns1[,i]~dt, col=(i+1))}



# Residual Autocorrelation
par(mfrow=c(2,length(models))); acfs <- list(); pacfs <- list()
for(i in 1:length(resids)){
  acfs[[i]] <- acf(resids[[i]], lag.max = (24*365.24), main=paste0('Residual Autocorrelation: ', modelNames[i]))
}
for(i in 1:length(resids)){
  pacfs[[i]] <- pacf(resids[[i]], lag.max = (24*4), main=paste0('Residual Partial Autocorrelation: ', modelNames[i]), ylim=c(-.1,0.2))
}

N=dim(dat)[1]
resMA24CubRegSeas=rollapply(resCubRegSeas, width = 24, FUN = mean, by = 24, fill = NA, align = "right")
resMA24CubRegSeas <- c(rep(0,24),rep(na.omit(resMA24CubRegSeas),each=24))[1:N]
updatedPredsCubRegSeas <- preds[[5]]+resMA24CubRegSeas
updatedResCubRegSeas <- dat$resNov-updatedPredsCubRegSeas

head(cbind(dat$resNov,
           preds[[5]],
           resCubRegSeas,
           resMA24CubRegSeas, 
           updatedPredsCubRegSeas,
           updatedResCubRegSeas
          ),100)

par(mfrow=c(2,3))
plot(dat$resNov~dt, type='l', lwd=3, col='darkgrey', main='PolyReg Prediction')
lines(preds[[5]]~dt, type='l', ylim=c(-.2,.7), col='blue')
plot(dat$resNov~dt, type='l', lwd=3, col='darkgrey', main='Corrected PolyReg Prediction')
lines(updatedPredsCubRegSeas~dt, type='l', col='blue', ylim=c(-.2,.7))
plot(dat$resNov~dt, type='l', lwd=3, col='darkgrey')
#lines(rollapply(resCubRegSeas, width = 24, FUN = mean, by = 1, fill = 0, align = "right")~dt, type='l', main='24hr rolling average', col='blue', ylim=c(-.2,.7), ylab='Rolling Average Prediction')
plot(resCubRegSeas~dt, type='l', ylim=c(-.2,.2), main='PolyReg Residuals')
plot(updatedResCubRegSeas~dt, type='l', ylim=c(-.2,.2), main='Corrected PolyReg Residuals')
plot(resMA24CubRegSeas~dt, type='l', ylim=c(-.2,.2), main='PolyReg Residuals Correction')

c(rep(0,24),rep(na.omit(resMA24CubRegSeas),each=24))[1:N]


residsMA12=list(resMA12LinReg=rollapply(resLinReg, width = 12, FUN = mean, by = 1, fill = 0, align = "right"),
                resMA12QuadReg=rollapply(resQuadReg, width = 12, FUN = mean, by = 1, fill = 0, align = "right"),
                resMA12CubReg=rollapply(resCubReg, width = 12, FUN = mean, by = 1, fill = 0, align = "right"),
                resvLinRegSeas=rollapply(resLinRegSeas, width = 12, FUN = mean, by = 1, fill = 0, align = "right"),
                resMA12CubRegSeas=rollapply(resCubRegSeas, width = 12, FUN = mean, by = 1, fill = 0, align = "right"))
residsMA7d=list(resMA7dLinReg=rollapply(resLinReg, width = (7*24), FUN = mean, by = 1, fill = 0, align = "right"),
                resMA7dQuadReg=rollapply(resQuadReg, width = (7*24), FUN = mean, by = 1, fill = 0, align = "right"),
                resMA7dCubReg=rollapply(resCubReg, width = (7*24), FUN = mean, by = 1, fill = 0, align = "right"),
                resMA7dLinRegSeas=rollapply(resLinRegSeas, width = (7*24), FUN = mean, by = 1, fill = 0, align = "right"),
                resMA7dCubRegSeas=rollapply(resCubRegSeas, width = (7*24), FUN = mean, by = 1, fill = 0, align = "right"))

# Residual Rolling Average
par(mfrow=c(1,length(models)))
for(i in 1:length(models)){
  plot(resids[[i]][1:2000]~dt[1:2000], type='l', ylim=c(-0.35, 0.2), 
       main=paste0(modelNames[i],'-',round(SSEs[i],2),'-(',minRes[i],',',maxRes[i],')'), 
       xlab='date', ylab='residuals',lwd=3)
  lines(residsMA12[[i]][12:2000]~dt[12:2000], lty=2, col='green',lwd=2)
  lines(residsMA7d[[i]][168:2000]~dt[168:2000], lty=2, col='red',lwd=2)
}

N <- dim(dat)[1]; updatedPreds <- list()
for(i in 1:length(models)){
  updatedPreds[[i]] <- rep(NA,N)
  updatedPreds[[i]][12:N] <- preds[[i]][12:N]+residsMA12[[i]][12:N]
}

for(i in 1:length(models)){
  plot(dat$resNov[1:480]~dt[1:480], type='l', lwd=2, ylim=c(-.2,0.65), main=modelNames[i], ylab='resNov', xlab='date')
  points(preds[[i]][1:480]~dt[1:480], type='p', col='blue',)
  points(updatedPreds[[i]][1:480]~dt[1:480], type='p', col='red',)
}

updatedResids <- list()
par(mfrow=c(1,length(resids)))
for(i in 1:(length(resids))){
  updatedResids[[i]] <- dat$resNov-updatedPreds[[i]]
  SSEs[i+length(resids)] <- round(sum(na.omit((updatedResids[[i]])^2)),4)
  MSEs[i+length(resids)] <- round(mean(na.omit((updatedResids[[i]])^2)),4)
  RMSEs[i+length(resids)] <- round(sqrt(mean(na.omit((updatedResids[[i]])^2))),4)
  minRes[i+length(resids)] <- round(min(na.omit((updatedResids[[i]]))),2)
  maxRes[i+length(resids)] <- round(max(na.omit((updatedResids[[i]]))),2)
  plot((updatedResids[[i]])[12:2000]~dt[12:2000], type='l', ylim=c(-0.4, 0.2), main=paste0(modelNames[i],'-',round(SSEs[i+length(resids)],2),'-(',minRes[i+length(resids)],',',maxRes[i+length(resids)],')'), xlab='date', ylab='updated residuals')
  abline(h=0, lty=2, col='grey')
}

par(mfrow=c(2,ceiling(length(models)/2)))
for(i in 1:length(models)){
  qqnorm(updatedResids[[i]], main=modelNames[i])
  qqline(updatedResids[[i]])
}

plot(residsMA12[[i]][12:2000]~dt[12:2000], lty=2, col='green')
plot(residsMA12[[i]][1:2000]~dt[1:2000], lty=2, col='regreend')


plot(residsMA7d[[i]][1:2000]~dt[1:2000], lty=2, col='red')
plot(residsMA7d[[i]][168:2000]~dt[168:2000], lty=2, col='red')

plot(residsMA7d[[i]][1:2000], lty=2, col='red')
