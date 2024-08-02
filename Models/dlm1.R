rm(list=ls())
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/Coastal Flood Predictions/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/Coastal Flood Predictions/Data/data_extract")

#Libraries
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(forecast)
library(dlm)
library(vars)
library(urca)
library(tseries)

X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)

X$datetime <- as.POSIXct(X$datetime, tz="GMT")
y_nov$datetime <- as.POSIXct(y_nov$datetime, tz="GMT")
y_pet$datetime <- as.POSIXct(y_pet$datetime, tz="GMT")
y_row$datetime <- as.POSIXct(y_row$datetime, tz="GMT")

head(X)
# Ocean_onshorewind: Ocean Wind, located at the NDBC buoy #46026,  On-shore component (from 100 N) 
# Gnoss_onshorewind: Local Wind at the Gnoss Field Airport, On-shore component (from 60 N) 
# AtmPres: Atmospheric Pressure, mBAR
# napa_flow_cfs: River flow at Napa Ricer at USGS gauge 11458000
head(y_nov)
head(y_pet)
head(y_row)
# Field residual =  stage - predicted It represents the variation of the water level due to non-tidal forcing.
# Field stage_m: raw water level data - detrended by removing the mean value, and resampled to hourly time intervals. Small data gaps were filled by linear interpolation. 
# Field predicted_m:  This is the predicted tide. The predicted tide was calculated using a publicly available Python routine based on a well documented Matlab routine called Utide (http://www.po.gso.uri.edu/~codiga/utide/utide.htm). 

dat_nov <- merge(X, y_nov, by='datetime')
dat_nov <- dat_nov[,c(1,2,3,4,5,8,7)]

#dat_nov <- dat_nov[1:1000,]
dt <- dat_nov[,1]
dat <- ts(data=dat_nov[2:7], start=c(2019, 20),
          end=c(2019,1019), frequency=365.25*24, class="matrix") #end=c(2021,12,28,15) #c(2022,6477)
dat <- as.data.frame(dat)
X <- dat[,1:5]

plot(dt, dat$stage_m, type='l')

Vinit <- var(dat$stage_m)

### DLM Regression + Seasonality
buildDLM <- function(params){
  s=12
  dlm1 <- dlmModReg(X) + dlmModSeas(s)
  
  V(dlm1) <- exp(params[1])
  diag(W(dlm1))[2:7] <- exp(params[2:7])
  
  return(dlm1)
}
fit <- dlmMLE(dat$stage_m, parm =c(log(Vinit), rep(0.5, ncol(dat)), 1), build = buildDLM)
dlm1 <- buildDLM(fit$par)
fit
dlm1F <- dlmFilter(dat$stage_m, dlm1); head(dlm1F$m,10)
# trend + 5 variables + s-1 seasonal

plot(dt, dat$stage_m, type='l')
lines(dt, dropFirst(dlm1F$m), col='red', lty=2)
lines(dt, dropFirst(dlm1F$m[,3]), col='red', lty=2)

par(mfrow=c(4,3))
plot(dt[150:1000], dat$stage_m[150:1000], type='l')
plot(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l')
for(i in 1:9){
  plot(dt[150:1000], dropFirst(dlm1F$m[,i])[150:1000], type='l', main=paste0(i))
}

### DLM Regression + Trig seasonality

buildDLM <- function(params){
  dlm1 <- dlmModReg(X) + dlmModTrig(s=(365.25*24), q=1)
  
  V(dlm1) <- exp(params[1])
  diag(W(dlm1))[2:8] <- exp(params[2:8])
  
  return(dlm1)
}
fit <- dlmMLE(dat$stage_m, parm =c(log(Vinit), rep(0.5, ncol(dat)+2)), build = buildDLM)
fit
dlm1 <- buildDLM(fit$par)
dlm1F <- dlmFilter(dat$stage_m, dlm1)

plot(dt, dat$stage_m, type='l')
lines(dt, dropFirst(dlm1F$m), col='red', lty=2)

par(mfrow=c(3,3))
plot(dt[150:1000], dat$stage_m[150:1000], type='l')
#plot(dt[1:1000], dropFirst(rowSums(dlm1F$m))[1:1000], type='l')
for(i in 1:8){
  plot(dt[150:1000], dropFirst(dlm1F$m[,i])[150:1000], type='l', main=paste0(i))
}
par(mfrow=c(1,1))
plot(dt[150:1000], dat$stage_m[150:1000], type='l', lwd=2, ylim=c(min(dat$stage_m),30))
lines(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l', col='red', lwd=2, lty=2)











#### RESIDUAL LEVEL

rm(list=ls())
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/Coastal Flood Predictions/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/Coastal Flood Predictions/Data/data_extract")

#Libraries
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(forecast)
library(dlm)
library(vars)
library(urca)
library(tseries)

X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)

X$datetime <- as.POSIXct(X$datetime, tz="GMT")
y_nov$datetime <- as.POSIXct(y_nov$datetime, tz="GMT")
y_pet$datetime <- as.POSIXct(y_pet$datetime, tz="GMT")
y_row$datetime <- as.POSIXct(y_row$datetime, tz="GMT")

head(X)
# Ocean_onshorewind: Ocean Wind, located at the NDBC buoy #46026,  On-shore component (from 100 N) 
# Gnoss_onshorewind: Local Wind at the Gnoss Field Airport, On-shore component (from 60 N) 
# AtmPres: Atmospheric Pressure, mBAR
# napa_flow_cfs: River flow at Napa Ricer at USGS gauge 11458000
head(y_nov)
head(y_pet)
head(y_row)
# Field residual =  stage - predicted It represents the variation of the water level due to non-tidal forcing.
# Field stage_m: raw water level data - detrended by removing the mean value, and resampled to hourly time intervals. Small data gaps were filled by linear interpolation. 
# Field predicted_m:  This is the predicted tide. The predicted tide was calculated using a publicly available Python routine based on a well documented Matlab routine called Utide (http://www.po.gso.uri.edu/~codiga/utide/utide.htm). 

dat_nov <- merge(X, y_nov, by='datetime')
dat_nov <- dat_nov[,c(1,2,3,4,5,6)]

dat_nov <- dat_nov[1:1000,]
dt <- dat_nov[,1]
dat <- ts(data=dat_nov[2:6], start=c(2019, 20),
          end=c(2019,1019), frequency=365.25*24, class="matrix") #end=c(2021,12,28,15) #c(2022,6477)
dat <- as.data.frame(dat)
X <- dat[,1:4]



### DLM Regression + Trig seasonality

buildDLM <- function(params){
  dlm1 <- dlmModReg(X) + dlmModTrig(s=(365.25*24), q=1)
  
  V(dlm1) <- exp(params[1])
  diag(W(dlm1))[2:7] <- exp(params[2:7])
  
  return(dlm1)
}
fit <- dlmMLE(dat$residual_m, parm =c(1, rep(0.5, ncol(dat)+2)), build = buildDLM)
fit
dlm1 <- buildDLM(fit$par)
dlm1F <- dlmFilter(dat$residual_m, dlm1)

plot(dt, dat$residual_m, type='l')
lines(dt, dropFirst(dlm1F$m), col='red', lty=2)

par(mfrow=c(3,3))
plot(dt[150:1000], dat$residual_m[150:1000], type='l')
plot(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l')
for(i in 1:7){
  plot(dt[150:1000], dropFirst(dlm1F$m[,i])[150:1000], type='l', main=paste0(i))
}
par(mfrow=c(1,1))
plot(dt[150:1000], dat$residual_m[150:1000], type='l', lwd=2, ylim=c(min(c(dat$residual_m, dropFirst(rowSums(dlm1F$m)))),max(c(dat$residual_m, dropFirst(rowSums(dlm1F$m))))))
lines(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l', col='red', lwd=2, lty=2)


### Simple Poly2 DLM with no regression or seasonality
buildDLM <- function(params){
  dlm1 <- dlmModPoly(2)
  
  V(dlm1) <- exp(params[1])
  diag(W(dlm1))[2] <- exp(params[2])
  
  return(dlm1)
}
fit <- dlmMLE(dat$residual_m, parm =c(1, rep(0.5, 1)), build = buildDLM)
fit
dlm1 <- buildDLM(fit$par)
dlm1F <- dlmFilter(dat$residual_m, dlm1)

plot(dt, dat$residual_m, type='l')
lines(dt, dropFirst(dlm1F$m), col='red', lty=2)

par(mfrow=c(2,2))
plot(dt[150:1000], dat$residual_m[150:1000], type='l')
plot(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l')
for(i in 1:2){
  plot(dt[150:1000], dropFirst(dlm1F$m[,i])[150:1000], type='l', main=paste0(i))
}
par(mfrow=c(1,1))
plot(dt[150:1000], dat$residual_m[150:1000], type='l', lwd=2, ylim=c(min(c(dat$residual_m, dropFirst(rowSums(dlm1F$m)))),max(c(dat$residual_m, dropFirst(rowSums(dlm1F$m))))))
lines(dt[150:1000], dropFirst(rowSums(dlm1F$m))[150:1000], type='l', col='red', lwd=2, lty=2)
