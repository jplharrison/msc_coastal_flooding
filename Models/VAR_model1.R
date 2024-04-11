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

diff_data <- dat_nov[2:dim(dat_nov)[1],]
diff_data$atmPres <- diff(dat_nov$atmPres)
diff_data$gnossWind <- diff(dat_nov$gnossWind)
diff_data$napaFlow <- diff(dat_nov$napaFlow)
diff_data$oceanWind <- diff(dat_nov$oceanWind)
diff_data$resNov <- diff(dat_nov$resNov)

par(mfrow=c(2,1))
plot(dat_nov$atmPres, type='l')
plot(diff_data$atmPres, type='l')

diff_ts <- ts(data=diff_data[,2:6], start=c(2019, 1, 1, 21),
              end=c(2022,9,27,15), deltat=(1/(365.25*24)), class="matrix")

train_window=1:8760; test_window=8761:9000;
train_data <- as.data.frame(scale(dat_nov[train_window,2:6], center = T, scale = T)); 
train_dt <- dat_nov[train_window,1]
train_dat <- ts(data=train_data[,1:5],  start=c(2019, 1, 1, 20),
                end=c(2020,1,1,19), frequency=365*24-1, class="matrix") #end=c(2021,12,28,15)

test_data <- as.data.frame(scale(dat_nov[test_window,2:6], center = T, scale = T))
test_dt <- dat_nov[test_window,1]
test_dat <-  ts(data=test_data[,1:4], start=c(2020, 1, 1, 20), end=c(2021,1,1,19), frequency=365*24-1, class="matrix")
test_resp <- ts(data=test_data[,5],   start=c(2020, 1, 1, 20), end=c(2021,1,1,19), frequency=365*24-1, class="matrix") 
test_dat <-  test_dat[1:240,]; #test_dat <- as.data.frame(test_dat)
test_resp <- test_resp[1:240]

adf.test(train_dat[,1])  #p-value = 0.01 - alternative hypothesis: stationary
adf.test(train_dat[,2])
adf.test(train_dat[,3])
adf.test(train_dat[,4])``
adf.test(train_dat[,5])

ggplot(data=cbind(train_dt, train_data)) + geom_point(mapping=aes(x=train_dt, y=resNov))

plot1 <- ggplot(data=train_data) + geom_point(mapping=aes(x=atmPres, y=resNov))
plot2 <- ggplot(data=train_data) + geom_point(mapping=aes(x=gnossWind, y=resNov))
plot3 <- ggplot(data=train_data) + geom_point(mapping=aes(x=napaFlow, y=resNov))
plot4 <- ggplot(data=train_data) + geom_point(mapping=aes(x=oceanWind, y=resNov))
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)

ind=which(train_data$napaFlow<1000)
ggplot(data=train_data[ind,]) + geom_point(mapping=aes(x=napaFlow, y=resNov))
ggplot(data=train_data[-ind,]) + geom_point(mapping=aes(x=napaFlow, y=resNov))


ols1 <- lm(train_resNov~train_atmPres+train_gnossWind+train_napaFlow+train_oceanWind)

par(mfrow=c(2,3))
acf(train_resNov, lag.max = 240)
acf(train_atmPres, lag.max = 8000)
acf(train_gnossWind, lag.max = 8000)
acf(train_napaFlow, lag.max = 8000)
acf(train_oceanWind, lag.max = 8000)

par(mfrow=c(2,3))
pacf(train_resNov, lag.max = 8000)
pacf(train_atmPres, lag.max = 8000)
pacf(train_gnossWind, lag.max = 8000)
pacf(train_napaFlow, lag.max = 8000)
pacf(train_oceanWind, lag.max = 8000)

# Finding optimal lags
var_select <- VARselect(train_dat, lag.max = 72, type = "const")
var_select$selection

# Building model
var1 <- VAR(train_dat, p=64, type='const', season=NULL, exog=NULL)
summ <- capture.output(summary(var1))
writeLines(summ, "var_summary.txt")

# Predictions 
train_pred <- predict(var1, n.ahead = 240)

par(mfrow=c(3,2))
plot(train_pred$fcst$resNov[,1], type='l', col=2, ylim=c(-0.04,0.12), lty=2, lwd=2, main='Novato Water Level')
lines(dat_nov[8761:9000,6], lwd=2)
plot(train_pred$fcst$atmPres[,1], type='l', col=2, lty=2, lwd=2, main='AtmPres', ylim=c(1017,1035))
lines(dat_nov[8761:9000,2], lwd=2)
plot(train_pred$fcst$gnossWind[,1], type='l', col=2, lty=2, lwd=2, main='Local Wind', ylim=c(-6,4))
lines(dat_nov[8761:9000,3], lwd=2)
plot(train_pred$fcst$napaFlow[,1], type='l', col=2, lty=2, lwd=2, main='River Flow')
lines(dat_nov[8761:9000,4], lwd=2)
plot(train_pred$fcst$oceanWind[,1], type='l', col=2, lty=2, lwd=2, main='Ocean Wind', ylim=c(-8,8))
lines(dat_nov[8761:9000,5], lwd=2)



lag_dat <- cbind(dat_nov$atmPres[1:(8760-12)], 
                 dat_nov$atmPres[2:(8760-11)], 
                 dat_nov$atmPres[3:(8760-10)], 
                 dat_nov$atmPres[4:(8760-9)], 
                 dat_nov$atmPres[5:(8760-8)], 
                 dat_nov$atmPres[6:(8760-7)], 
                 dat_nov$atmPres[7:(8760-6)], 
                 dat_nov$atmPres[8:(8760-5)], 
                 dat_nov$atmPres[9:(8760-4)], 
                 dat_nov$atmPres[10:(8760-3)], 
                 dat_nov$atmPres[11:(8760-2)],  
                 dat_nov$atmPres[12:(8760-1)], 
                 dat_nov$atmPres[13:(8760)],
                 dat_nov$gnossWind[1:(8760-12)], 
                 dat_nov$gnossWind[2:(8760-11)], 
                 dat_nov$gnossWind[3:(8760-10)], 
                 dat_nov$gnossWind[4:(8760-9)], 
                 dat_nov$gnossWind[5:(8760-8)], 
                 dat_nov$gnossWind[6:(8760-7)], 
                 dat_nov$gnossWind[7:(8760-6)], 
                 dat_nov$gnossWind[8:(8760-5)], 
                 dat_nov$gnossWind[9:(8760-4)], 
                 dat_nov$gnossWind[10:(8760-3)], 
                 dat_nov$gnossWind[11:(8760-2)],  
                 dat_nov$gnossWind[12:(8760-1)], 
                 dat_nov$gnossWind[13:(8760)],
                 dat_nov$napaFlow[1:(8760-12)], 
                 dat_nov$napaFlow[2:(8760-11)], 
                 dat_nov$napaFlow[3:(8760-10)], 
                 dat_nov$napaFlow[4:(8760-9)], 
                 dat_nov$napaFlow[5:(8760-8)], 
                 dat_nov$napaFlow[6:(8760-7)], 
                 dat_nov$napaFlow[7:(8760-6)], 
                 dat_nov$napaFlow[8:(8760-5)], 
                 dat_nov$napaFlow[9:(8760-4)], 
                 dat_nov$napaFlow[10:(8760-3)], 
                 dat_nov$napaFlow[11:(8760-2)],  
                 dat_nov$napaFlow[12:(8760-1)], 
                 dat_nov$napaFlow[13:(8760)],
                 dat_nov$oceanWind[1:(8760-12)], 
                 dat_nov$oceanWind[2:(8760-11)], 
                 dat_nov$oceanWind[3:(8760-10)], 
                 dat_nov$oceanWind[4:(8760-9)], 
                 dat_nov$oceanWind[5:(8760-8)], 
                 dat_nov$oceanWind[6:(8760-7)], 
                 dat_nov$oceanWind[7:(8760-6)], 
                 dat_nov$oceanWind[8:(8760-5)], 
                 dat_nov$oceanWind[9:(8760-4)], 
                 dat_nov$oceanWind[10:(8760-3)], 
                 dat_nov$oceanWind[11:(8760-2)],  
                 dat_nov$oceanWind[12:(8760-1)], 
                 dat_nov$oceanWind[13:(8760)],
                 dat_nov$resNov[1:(8760-12)], 
                 dat_nov$resNov[2:(8760-11)], 
                 dat_nov$resNov[3:(8760-10)], 
                 dat_nov$resNov[4:(8760-9)], 
                 dat_nov$resNov[5:(8760-8)], 
                 dat_nov$resNov[6:(8760-7)], 
                 dat_nov$resNov[7:(8760-6)], 
                 dat_nov$resNov[8:(8760-5)], 
                 dat_nov$resNov[9:(8760-4)], 
                 dat_nov$resNov[10:(8760-3)], 
                 dat_nov$resNov[11:(8760-2)],  
                 dat_nov$resNov[12:(8760-1)], 
                 dat_nov$resNov[13:(8760)]
)
variables <- colnames(dat_nov)[2:6]
names <- c()
for(i in 1:5){
  for(j in 0:11){
    names <- c(names, paste0(variables[i],'-',12-j))
  }
  names <- c(names, paste0(variables[i]))
}
colnames(lag_dat) <- names
#lag_dat <- as.data.frame(lag_dat)

lagOLS <- lm(resNov~., data=lag_dat)


ols_pred <- predict(ols1, newdata = test_dat)
ols_pred <- ols1$coefficients[1]+ols1$coefficients[2]*test_dat[,1]+
  ols1$coefficients[3]*test_dat[,2]+ols1$coefficients[4]*test_dat[,3]+
  ols1$coefficients[5]*test_dat[,4]
ols_predIS <- ols1$coefficients[1]+ols1$coefficients[2]*train_dat[,1]+
  ols1$coefficients[3]*train_dat[,2]+ols1$coefficients[4]*train_dat[,3]+
  ols1$coefficients[5]*train_dat[,4]

par(mfrow=c(2,1))
plot(test_resp, type='l', lwd=2, ylim=c(min(c(test_resp,ols_pred)),max(c(test_resp,ols_pred))))
lines(ols_pred, type='l', lwd=2, lty=2, col=2)
plot(train_dat[,5], type='l', lwd=2, ylim=c(min(c(train_dat[,5],ols_predIS)),max(c(train_dat[,5],ols_predIS))))
lines(ols_predIS, type='l', lwd=2, lty=2, col=2)


lag_pred <- lag_dat






