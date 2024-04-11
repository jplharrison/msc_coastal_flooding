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

# Merge predictors and variables into data frame
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
dat <- as.data.frame(dat)

# OLS Multiple Linear Regression model
ols <- lm(resNov~., data = dat)
# Fitted Values
resNovIS <- predict(ols, newdata=dat[,1:4])
# Insample MSE
MSEIS  <- sum((resNovIS - dat[,5])^2)
# Insample Residual Errors
resIS <- dat[,5] - resNovIS
resIS <- ts(data=resIS, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

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


