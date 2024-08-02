rm(list=ls())

# Data Directory
#setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")

#Libraries
library(ggplot2)
library(gridExtra)
library(forecast)
library(vars)
library(urca)
library(tseries)
library(tidyverse)

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
tm_14 <- time(dat); start_tm <- min(tm_14); end_tm <- max(tm_14); yrs <- (end_tm-start_tm)*365.25/14; tm_14 <- ((tm_14-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
dat <- as.data.frame(cbind(sin(tm_14),cos(tm_14),sin(tm_28),cos(tm_28),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:6] <- c('sin_14','cos_14','sin_28','cos_28','sin_yr','cos_yr')

# OLS Multiple Linear Regression model
ols <- lm(resNov~atmPres+gnossWind+napaFlow+oceanWind, data = dat)
# Fitted Values
resNovIS <- predict(ols, newdata=dat[,7:10])
# Insample MSE
MSEIS  <- sum((resNovIS - dat$resNov)^2)
# Insample Residual Errors
resIS <- ts(data=resid(ols) , start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

anova(ols)


ols2 <-  lm(resNov~0+sin_14+cos_14+sin_28+cos_28+sin_yr+cos_yr+
             atmPres+gnossWind+napaFlow+oceanWind
           , data = dat); summary(ols2)
#resNovIS2 <- predict(ols2, newdata=dat[,5:9])
resNovIS2 <- predict(ols2, newdata=dat[,1:10])
resIS2 <- resid(ols2)
resIS2 <- ts(data=resIS2, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")
# Plot observed, fitted and residuals
par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$resNov, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(resNovIS2, dat$resNov)),max(c(resNovIS2, dat$resNov))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(resNovIS2, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS2, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS2, resIS2)),max(c(resIS2, resIS2))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')
par(mfrow=c(1,2))
qqnorm(resIS2); qqline(resIS2)
hist(resIS2, breaks=30, main='Distribution of Residuals/Errors',probability = F)
par(mfrow=c(1,1))
a=acf(resIS2, lag.max = 8760, plot = T, main='ACF Residuals/Errors - ols2 w/ trig seas')
abline(v=0:26*(14*24)/(365.25*24), lty=2, col='grey')
abline(v=0:13*(27*24+8)/(365.25*24), lty=2, col='grey')
lines((1:8760)/8760, dat$cos_28[1:8760], type='l', lty=2, col='blue')
plot(dt, resIS2, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS2, resIS2)),max(c(resIS2, resIS2))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(dt, dat$cos_28, type='l', lty=2, col='blue')

par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$resNov, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(resNovIS2, dat$resNov)),max(c(resNovIS2, dat$resNov))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(resNovIS2, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$resNov, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(resNovIS, dat$resNov)),max(c(resNovIS, dat$resNov))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(resNovIS, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)


# Plot observed, fitted and residuals
par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat[,5], type='l',lwd=2, ylab='resNov Response', xlab='Date', xaxt='n', ylim=c(min(c(resNovIS, dat[,5])),max(c(resNovIS, dat[,5]))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(resNovIS, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS, type='l',lwd=2, ylab='resNov Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')

# Distribution of residual errors
par(mfrow=c(1,2))
qqnorm(resIS)
hist(resIS, breaks=30, main='Distribution of Residuals/Errors',probability = F)

adf.test(resIS)

# ACF of Residual errors
par(mfrow=c(1,1))
a=acf(resIS, lag.max = 8760, plot = T, main='ACF Residuals/Errors')
cutoff=0.05
abline(h=c((-1*cutoff),cutoff), lty=2, col='red', lwd=2)
abline(v=1:8*(100*24)/(365.25*24), lty=2, col='grey')
abline(v=240/(365.25*24), lty=2, col='grey')
points((which(abs(a$acf)>cutoff)/(365.25*24)), a$acf[which(abs(a$acf)>cutoff)], col='green', pch=19)

# ARIMA model for Residual Error Correction
arima1 <- auto.arima(resIS, seasonal = T)
arima1_forecast <- forecast(arima1, 240, level=c(50,95))
autoplot(arima1_forecast)

arima2 <- Arima(resIS, order=c(240, 0, 0), seasonal=c(0.5,0,0))
arima2_forecast <- forecast(arima2, h=240, level=c(50,95))
autoplot(arima1_forecast)


arima3 <- Arima(...)

arima3_forecast <- forecast(arima3, 240, level=c(50,95))
autoplot(arima3_forecast)

# Explore seasonal component of residual errors
dec <- decompose(resIS)
des <- dec$x-dec$seasonal

par(mfrow=c(3,2))
plot(dec$x, type='l', ylim=c(min(dec$x), max(dec$x))); 
plot(des, type='l')
plot(ma(dec$x, order=(24*7)), lty=2, col='grey', ylim=c(min(dec$x), max(dec$x))) #lines(dec$x-dec$seasonal, lty=2, col='grey')
plot(ma(des, order=(24*7)), lty=2, col='grey', ylim=c(min(dec$x), max(dec$x))) 
plot(ma(dec$x, order=(24*30)),type='b', lty=2, col='grey', ylim=c(min(dec$x), max(dec$x))) #lines(dec$x-dec$seasonal, lty=2, col='grey')
plot(ma(des, order=(24*30)), type='b', lty=2, col='grey', ylim=c(min(dec$x), max(dec$x))) 


