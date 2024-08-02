
rm(list=ls())

# Data Directory
#setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")

#Libraries
library(ggplot2)
library(gridExtra)
library(forecast)
library(vars)
library(urca)
library(tseries)
library(tidyverse)
library(splines)
library(visreg)
library(car)
library(mgcv)
library(MASS)

# Load ata
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)
# Standardise - centre and scale data - not necessary
#X[,2:5] <- scale(X[,2:5], center = T, scale=T)
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
# Field stage_m: raw water level data - detrended by removing the [overall] mean value, and resampled to hourly time intervals. Small data gaps were filled by linear interpolation. 
# Field predicted_m:  This is the predicted tide. The predicted tide was calculated using a publicly available Python routine based on a well documented Matlab routine called Utide (http://www.po.gso.uri.edu/~codiga/utide/utide.htm). 

#### Include predicted as an independent variable; use stage_m as the response.
dat_nov <- merge(X, y_nov, by='datetime')
dat_nov <- dat_nov[,c(1,2,3,4,5,8,7)]

#dat_nov <- dat_nov[1:1000,]
dt <- dat_nov[,1]
dat <- ts(data=dat_nov[2:7], start=c(2019, 20),
          end=c(2022,6477), frequency=365.25*24, class="matrix") #end=c(2021,12,28,15) #
tm_yr <- time(dat); start_tm <- min(tm_yr); end_tm <- max(tm_yr); yrs <- end_tm-start_tm; tm_yr <- ((tm_yr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
tm_28 <- time(dat); start_tm <- min(tm_28); end_tm <- max(tm_28); yrs <- (end_tm-start_tm)*365.25/27.33; tm_28 <- ((tm_28-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
tm_qr <- time(dat); start_tm <- min(tm_qr); end_tm <- max(tm_qr); yrs <- (end_tm-start_tm)*4; tm_qr <- ((tm_qr-start_tm)/(end_tm-start_tm)) * (yrs*2*pi)
#dat <- as.data.frame(cbind(tm_yr,tm_28,dat)); colnames(dat)[1:2] <- c('tm_yr','tm_28')
dat <- as.data.frame(cbind(sin(tm_28),cos(tm_28),sin(tm_qr),cos(tm_qr),sin(tm_yr),cos(tm_yr),dat)); colnames(dat)[1:6] <- c('sin_28','cos_28','sin_qr','cos_qr','sin_yr','cos_yr')
#new_dat <- as.data.frame(cbind(dat[,1:5],time(dat))); colnames(new_dat)[6] <- 'tm'



# OLS Multiple Linear Regression model
ols <- lm(stage_m~0+AtmPres+Gnoss_onshorewind+napa_flow_cfs+ocean_onshorewind+predicted #+bs(dt,) add cyclic spline for seasonal component 
          , data = dat); summary(ols)
ols1 <-  gam(stage_m~0+AtmPres+Gnoss_onshorewind+napa_flow_cfs+ocean_onshorewind+predicted + 
               s(tm_yr, bs='cc', k=10) +
               s(tm_28, bs='cc', k=10)
            , data = dat); summary(ols1)
main_lbl <- 'Regular OLS'; main_lbl1 <- 'GAM w/ 28day and 1year smooths'; main_lbl_test <- 'OLS w/ ns(predicted) + 28day + year seas';  

rrX <- dat[,1:9]
rry <- dat$stage_m
lmr <- lm.ridge(stage_m~0+sin_28+cos_28+sin_yr+cos_yr+
                  AtmPres+Gnoss_onshorewind+napa_flow_cfs+ocean_onshorewind+predicted
                , data = dat, lambda=100)
lmr.res <- ts(data=resid(lmr), start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

acf(lmr.res, lag.max = 8766)

ols_test1 <-  lm(stage_m~0+sin_28+cos_28+sin_qr+cos_qr+sin_yr+cos_yr+
                  AtmPres+Gnoss_onshorewind+napa_flow_cfs+ocean_onshorewind+predicted
                , data = dat); summary(ols_test1)
#stageNovIS_test1 <- predict(ols_test1, newdata=dat[,5:9])
resIS_test1 <- resid(ols_test1)
resIS_test1 <- ts(data=resIS_test1, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")


ols_test <-  lm(stage_m~0+sin_28+cos_28+sin_qr+cos_qr+sin_yr+cos_yr+
            AtmPres+Gnoss_onshorewind+napa_flow_cfs+ocean_onshorewind+ns(predicted, df = 4)
             , data = dat); summary(ols_test)
#stageNovIS_test <- predict(ols_test, newdata=dat[,5:9])
stageNovIS_test <- predict(ols_test, newdata=dat[,1:9])
resIS_test <- resid(ols_test)
resIS_test <- ts(data=resIS_test, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")
# Plot observed, fitted and residuals
par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$stage_m, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(stageNovIS_test, dat[,6])),max(c(stageNovIS_test, dat[,6]))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(stageNovIS_test, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS_test, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS_test, resIS_test)),max(c(resIS_test, resIS_test))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')
par(mfrow=c(1,2))
qqnorm(resIS_test); qqline(resIS_test)
hist(resIS_test, breaks=30, main='Distribution of Residuals/Errors',probability = F)
par(mfrow=c(1,1))
a=acf(resIS_test, lag.max = 8760, plot = T, main='ACF Residuals/Errors - ols_test w/ trig seas + pred spline')
abline(v=0:13*(27*24+8)/(365.25*24), lty=2, col='grey')
lines((1:8760)/8760, dat$cos_28[1:8760], type='l', lty=2, col='blue')
lines((1:8760)/8760, dat$cos2_28[1:8760], type='l', lty=2, col='red')
plot(dt, resIS_test, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS_test, resIS_test)),max(c(resIS_test, resIS_test))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
#lines(dt, dat$cos_28, type='l', lty=2, col='blue')
#lines((1:8760)/8760, dat$cos2_28[1:8760], type='l', lty=2, col='red')

par(mfrow=c(1,2))
plot(dt, resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS),max(resIS)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)
plot(dt, resIS_test, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS_test),max(resIS_test)), main=paste0('In Sample Residuals - ', main_lbl_test)); axis(1, at = ticks, labels = dt[ticks], las=1)

par(mfrow=c(1,2))
plot(dt, resIS-resIS_test, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS-resIS_test),max(resIS-resIS_test)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)
plot(dt, resIS-resIS_test1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS-resIS_test1),max(resIS-resIS_test1)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(v=dt[(1:52)*(27*24+8)], lty=2, col='grey')
abline(v=dt[(1:4)*(365.25*24)], lty=2, col='grey', lwd=3)

par(mfrow=c(1,2))
plot(dt, resIS-resIS_test1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS-resIS_test1),max(resIS-resIS_test1)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)
resIS_test1 <- as.matrix(dat[,1:11])%*%matrix(c(ols_test1$coefficients[1:6]*10, ols_test1$coefficients[7:11]), ncol=1)
plot(dt, resIS-resIS_test1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS-resIS_test1),max(resIS-resIS_test1)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)

par(mfrow=c(1,2))
plot(dt, resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS),max(resIS)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)
plot(dt, resIS_test1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS_test1),max(resIS_test1)), main=paste0('In Sample Residuals - ', main_lbl)); axis(1, at = ticks, labels = dt[ticks], las=1)


plot(ols1, select = 1, shade = TRUE, rug = TRUE,
     main = "Contribution of Cyclic Spline",
     xlab = "Time",
     ylab = "Effect on stage_m")
plot(ols1, select = 2, shade = TRUE, rug = TRUE,
     main = "Contribution of Cyclic Spline",
     xlab = "Time",
     ylab = "Effect on stage_m")

# Fitted Values
stageNovIS <- predict(ols, newdata=dat[,5:9])
stageNovIS1 <- predict(ols1, newdata=dat)
# Insample MSE
MSEIS  <- sum((stageNovIS - dat$stage_m)^2)
# Insample Residual Errors
resIS <- resid(ols)
resIS <- ts(data=resIS, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")
resIS1 <- resid(ols1)
resIS1 <- ts(data=resIS1, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")

# Plot observed, fitted and residuals
par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$stage_m, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(stageNovIS, dat[,6])),max(c(stageNovIS, dat[,6]))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(stageNovIS, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')

par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$stage_m, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(stageNovIS1, dat[,6])),max(c(stageNovIS1, dat[,6]))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(stageNovIS1, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS1),max(resIS1)), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')

par(mfrow=c(1,2))
plot(resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')
plot(resIS1, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(resIS1),max(resIS1)), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')

par(mfrow=c(1,1))
plot(dat$stage_m, type='l')#, ylim=c(-10,8))
#abline(h=ols$coefficients[1], col=2)
for(i in 1:9){
  print(colnames(dat)[i])
  lines(dat[,i]*ols$coefficients[i], col=i+2)
  Sys.sleep(5)
}


par(mfrow=c(2,2))
qqnorm(resIS); qqline(resIS)
hist(resIS, breaks=30, main='Distribution of Residuals/Errors',probability = F)

qqnorm(resIS1); qqline(resIS1)
hist(resIS1, breaks=30, main='Distribution of Residuals/Errors',probability = F)

### Using ns(predicted) - increases effect of AtmPres, slight decrease in effect of LocalWind
### effect of predicted increases with value of predicted
## note: effect := magnitude of coefficient
ols2 <- lm(stage_m~0+AtmPres+Gnoss_onshorewind+napa_flow_cfs+
             ocean_onshorewind+ns(predicted, df = 4), data = dat); summary(ols2)
stageNovIS2 <- predict(ols2, newdata=dat[,1:5])
# Insample Residual Errors
resIS2 <- resid(ols2)
resIS2 <- ts(data=resIS2, start=c(2019, 20), end=c(2022,6477), frequency=365.25*24, class="matrix")
par(mfrow=c(2,1))
ticks=floor(seq(from=1, to=dim(dat)[1], length=7))[c(2,4,6)]
plot(dat$stage_m, type='l',lwd=2, ylab='Nov stage_m', xlab='Date', xaxt='n', ylim=c(min(c(stageNovIS2, dat[,6])),max(c(stageNovIS2, dat[,6]))), main='In Sample Prediction'); axis(1, at = ticks, labels = dt[ticks], las=1)
lines(stageNovIS2, type='l', lty=2, lwd=2, col=2)
legend('topright', legend=c('Observed', 'Predicted'), bty='n', lty=c(1,2), col=c(1,2),  lwd=2)
plot(resIS2, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS2, resIS2)),max(c(resIS2, resIS2))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')

par(mfrow=c(1,2))
qqnorm(resIS2); qqline(resIS2);
hist(resIS2, breaks=30, main='Distribution of Residuals/Errors (OLS w/ Spline)',probability = F)


ols3 <- lm(stage_m~0+AtmPres+Gnoss_onshorewind+napa_flow_cfs+
             ocean_onshorewind+ns(predicted, df = 4)+
             napa_flow_cfs*ns(predicted, df = 4), data = dat); summary(ols3)
stageNovIS3 <- predict(ols3, newdata=dat[,1:5])
# Insample Residual Errors
resIS3 <- resid(ols3)

par(mfrow=c(1,2))
plot(resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=(-3:3)*0.2, lty=2, col='lightgrey')
plot(resIS2, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=(-3:3)*0.2, lty=2, col='lightgrey')
plot(resIS3, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=(-3:3)*0.2, lty=2, col='lightgrey')


# ANOVAs for original OLS model, OLS w/ Spline, OLS w/ Spline + NapaFlow:Predicted interaction
anova(ols)
anova(ols2)
anova(ols3)

par(mfrow=c(1,1))
plot(resIS, type='l',lwd=2, ylab='Nov stage_m Residuals/Errors', xlab='Date', xaxt='n', ylim=c(min(c(resIS, resIS)),max(c(resIS, resIS))), main='In Sample Residuals'); axis(1, at = ticks, labels = dt[ticks], las=1)
abline(h=0, lty=2, col='lightgrey')
lines(resIS-resIS2, col=2)
lines(resIS-resIS3, col=3)

par(mfrow=c(1,1))
plot(resIS2-resIS3, type='l')



### Scatterplots of stage_m and predictors with line of best fit
par(mfrow=c(2,2))
m <- lm(stage_m~AtmPres, data=dat)
plot(dat$AtmPres, dat$stage_m, pch=19)
abline(m$coefficients[1],m$coefficients[2], lty=2, col='blue', lwd=2)

m <- lm(stage_m~Gnoss_onshorewind , data=dat)
plot(dat$Gnoss_onshorewind , dat$stage_m, pch=19)
abline(m$coefficients[1],m$coefficients[2], lty=2, col='blue', lwd=2)

m <- lm(stage_m~napa_flow_cfs  , data=dat)
plot(dat$napa_flow_cfs  , dat$stage_m, pch=19)
abline(m$coefficients[1],m$coefficients[2], lty=2, col='blue', lwd=2)

m <- lm(stage_m~ocean_onshorewind  , data=dat)
plot(dat$ocean_onshorewind  , dat$stage_m, pch=19)
abline(m$coefficients[1],m$coefficients[2], lty=2, col='blue', lwd=2)


# 
par(mfrow=c(1,2))
adf.test(resIS)
adf.test(resIS2)

# ACF of Residual errors
par(mfrow=c(1,2))
a=acf(resIS, lag.max = 8760, plot = T, main='ACF Residuals/Errors - OLS')
abline(v=0:13*(27*24+8)/(365.25*24), lty=2, col='grey')# ACF w/ spline
#par(mfrow=c(1,1))
a=acf(resIS1, lag.max = 8760, plot = T, main='ACF Residuals/Errors - OLS w/ cyclic spline')
abline(v=0:13*(27*24+8)/(365.25*24), lty=2, col='grey')
par(mfrow=c(1,1))
a=acf(resIS2, lag.max = 8760, plot = T, main='ACF Residuals/Errors - OLS w/ spline')
abline(v=0:13*(27*24+8)/(365.25*24), lty=2, col='grey')

pacf(resIS, lag.max = 8760)
# Lunar cycle seasonality + smaller peaks inbetween 

# ARIMA model for Residual Error Correction
arima1 <- auto.arima(resIS, seasonal = T)
arima1_forecast <- forecast(arima1, 240, level=c(50,95))
autoplot(arima1_forecast)

arima2 <- Arima(resIS, order=c(240, 0, 0), seasonal=c(27*24,0,0))
arima2_forecast <- forecast(arima2, h=240, level=c(50,95))
autoplot(arima1_forecast)







