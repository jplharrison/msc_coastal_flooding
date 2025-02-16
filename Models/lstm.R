rm(list=ls())
#rm(X)
# Data Directory
setwd("C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")
#setwd("C:/Users/User/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Data/data_extract")

#Libraries
library(forecast)
library(ggplot2)
library(gridExtra)
library(nlme)
library(scales)
library(tidyverse)
library(tseries)
library(keras)
library(tensorflow)
library(tfruns)

# save.image(file='lstm_tuning.R')
# save.image(file='lstm24.R')
# save.image(file='lstm3.R')
# save.image(file='lstm3_rescale.R')
# save.image(file='lstm24_rescale.R')
# load(file='lstm24_rescale.R')
# load(file='lstm3_rescale.R')

# Load Data
X <- read.csv(file='Input_trainingset.csv', header=T)
y_nov <- read.csv(file='novato_wl_1hr_up.csv', header=T)
y_pet <- read.csv(file='petaluma_wl_1hr.csv', header=T)
y_row <- read.csv(file='rowland_wl_1hr.csv', header=T)

# Rescale predictors between 0 and 1
for(i in 2:5){X[,i] <- rescale(X[,i], to=c(-1,1))}
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

N <- dim(X)[1]; N_test <- 2448; N_train <- N-N_test
dat <- dat_nov %>% dplyr::select(atmPres, gnossWind, napaFlow, oceanWind, predNov, stageNov)

dt_train <- dat_nov$datetime[1:N_train]
dt_test <- dat_nov$datetime[(1+N_train):N]
dat_train <- ts(data=dat[1:N_train,], start=c(2019, 20),
                end=c(2022,(6477.562-N_test)), frequency=365.24219*24, class="matrix") 
dat_test <- as.data.frame(ts(data=dat[(1+N_train):N,], start=c(2022,(6477.562-N_test+1)),
                             end=c(2022,(6477.562)), frequency=365.24219*24, class="matrix"))
dat_train <- as.data.frame(dat_train); stage_test <- dat_test$stageNov; dat_test <- as.data.frame(dat_test[,-10]);

var_names <- c('Atmospheric Pressure', 'Local Wind', 'Napa Flow', 'Ocean Wind', 'Predicted Tide', 'Water Level')


# Prepare data with lags
prepare_data <- function(data, response, predictors, timesteps, forecast_horizon) {
  N_samples <- nrow(data)-timesteps-forecast_horizon+1
  
  x_data <- matrix(NA, nrow = N_samples, ncol = timesteps + length(predictors))
  y_data <- matrix(NA, nrow = N_samples, ncol = forecast_horizon)
  
  # Iterate over the rows to create lagged sequences
  for (i in 1:N_samples) {
    if(i%%500==0){print(i)}
    x_data[i,] <- unlist(c(data[i:(i+timesteps-1), response], data[i+timesteps,predictors]))
    y_data[i,] <- data[(i+timesteps):(i+timesteps+forecast_horizon-1), response]
  }
  
  #colnames(x_data) <- c(paste0('lag',24:1), predictors)
  
  x_array <- array(x_data, dim = c(nrow(x_data), 1, ncol(x_data)))
  y_array <- array(y_data, dim = c(nrow(y_data), forecast_horizon, 1))
  
  return(list(x = x_array, y = y_array))
}

# Specify parameters
predictors <- c("atmPres", "gnossWind", "napaFlow", "oceanWind", "predNov")
response <- "stageNov"
timesteps <- 3
forecast_horizon <- 96
# Prepare training and testing datasets 16:33:00
train_set <- prepare_data(dat_train, response, predictors, timesteps, forecast_horizon)
test_set <- prepare_data(dat_test, response, predictors, timesteps, forecast_horizon)

x_train <- train_set$x
y_train <- train_set$y
x_test <- test_set$x
y_test <- test_set$y

# colnames(x_train) <- c(paste0('lag',24:1), predictors)
# colnames(x_test) <- c(paste0('lag',24:1), predictors)
# colnames(y_train) <- response
# colnames(y_test) <- response

dim(x_train)
dim(y_train)

# model <- keras_model_sequential() %>%
#   layer_lstm(units = 50, input_shape = c(dim(x_train)[2:3]), return_sequences = TRUE) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_lstm(units = 50, return_sequences = FALSE) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = forecast_horizon)
# 
# # Compile the model
# model %>% compile(
#   loss = "mean_squared_error",
#   optimizer = "adam"
# )
# 
# # Train the model
# history <- model %>% fit(
#   x_train, y_train,
#   epochs = 50,
#   batch_size = 32,
#   validation_split = 0.2,
#   verbose = 1
# )
# 
# 
# # Make predictions
# predictions <- model %>% predict(x_test)
#
# resids <- stage_test[(1+timesteps):(2448-forecast_horizon+1)]-predictions
# 
# par(mfrow=c(3,4))
# for(i in 1:12){
#   pred_window <- (i+timesteps):(i+forecast_horizon+forecast_horizon-1)
#   plot(stage_test[pred_window]~dt_test[pred_window], type='l',lwd=3,col='lightgrey', main=paste0('Iter ',i))
#   lines(predictions[i,]~dt_test[pred_window], col='blue')
# }



### TUNING

# Define the hyperparameter grid
grid <- expand.grid(
  units = c(32, 50, 64),        # Number of LSTM units
  dropout = c(0.05, 0.15, 0.25),        # Dropout rate
  layers = c(1, 2, 3),             # Number of LSTM layers
  batch_size = c(32, 64),       # Batch size
  learning_rate = c(0.001) # Learning rate for optimizer
)

grid <- expand.grid(
  units = c(64, 96, 128),        # Number of LSTM units
  dropout = c(0, 0.05),        # Dropout rate
  layers = c(3,4),             # Number of LSTM layers
  batch_size = c(64),       # Batch size
  learning_rate = c(0.001) # Learning rate for optimizer
)

grid <- expand.grid(
  units = c(64),        # Number of LSTM units
  dropout = c(0.05),        # Dropout rate
  layers = c(2),             # Number of LSTM layers
  batch_size = c(32),       # Batch size
  learning_rate = c(0.001, 0.01, 0.1) # Learning rate for optimizer
)


build_model <- function(units, dropout, layers, input_shape, forecast_horizon, learning_rate) {
  model <- keras_model_sequential()
  
  # Add the first LSTM layer
  model %>%
    layer_lstm(units = units, input_shape = input_shape, return_sequences = (layers > 1)) %>%
    layer_dropout(rate = dropout)
  
  # Add additional LSTM layers if specified
  if (layers > 1) {
    model %>%
      layer_lstm(units = units, return_sequences = FALSE) %>%
      layer_dropout(rate = dropout)
  }
  
  # Add dense output layer
  model %>% layer_dense(units = forecast_horizon)
  
  # Compile the model
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adam(learning_rate = learning_rate)
  )
  
  return(model)
}

results <- list() # To store results
forecast_horizon <- dim(y_train)[2] # Output dimension
input_shape <- c(dim(x_train)[2], dim(x_train)[3]) # Input shape (timesteps, features)

for (i in seq_len(nrow(grid))) {
  params <- grid[i, ]
  
  cat("Training model with params: ", paste(names(params), params, sep=" = ", collapse = ", "), "\n")
  
  
  model <- build_model(
    units = params$units,
    dropout = params$dropout,
    layers = params$layers,
    input_shape = input_shape,
    forecast_horizon = forecast_horizon,
    learning_rate = params$learning_rate
  )
  
  history <- model %>% fit(
    x_train, y_train,
    epochs = 50,
    batch_size = params$batch_size,
    validation_split = 0.2,
    verbose = 1
  )
  
  # Evaluate on validation set
  val_loss <- min(history$metrics$val_loss)
  
  # Save results
  results[[i]] <- list(
    params = params,
    val_loss = val_loss
  )
}

# Find the best parameters
best_model <- results[[which.min(sapply(results, function(x) x$val_loss))]]
cat(paste0("Best hyperparameters: ", best_model$params, "\n"))

# results_save <- results
# best_model_save <- best_model
# best_params_save <- best_model$params

best_params <- best_model$params

final_model <- build_model(
  units = best_params$units,
  dropout = best_params$dropout,
  layers = best_params$layers,
  input_shape = input_shape,
  forecast_horizon = forecast_horizon,
  learning_rate = best_params$learning_rate
)

history <- final_model %>% fit(
  x_train, y_train,
  epochs = 50,
  batch_size = best_params$batch_size,
  validation_split = 0.2,
  verbose = 1
)

# Make predictions
train_preds <- final_model %>% predict(x_train)
predictions <- final_model %>% predict(x_test)

N_samples <- dim(x_train)[1]
train_stage <- matrix(NA, nrow=N_samples, ncol=forecast_horizon)
for(f in 1:forecast_horizon){train_stage[,f] <- dat_train$stageNov[(timesteps+f):(N_samples+timesteps+f-1)]}
train_resids <- t(train_stage-train_preds)
rmse <- function(x){return(sqrt(mean(x^2)))}
mae <- function(x){return(mean(abs(x)))}

mean(apply(train_resids, 2, rmse))

par(mfrow=c(3,4)); jump=48
for(i in (1+jump):(12+jump)){
  pred_window <- (i+timesteps):(i+timesteps+forecast_horizon-1)
  plot(train_stage[i,]~dt_test[pred_window], type='l',lwd=3,col='lightgrey', main=paste0('Iter ',i))
  lines(train_preds[i,]~dt_test[pred_window], col='blue')
}

test_stage <- matrix(NA, nrow=dim(predictions)[1], ncol=forecast_horizon)
for(f in 1:forecast_horizon){test_stage[,f] <- dat_test$stageNov[(timesteps+f):(dim(predictions)[1]+timesteps+f-1)]}
test_resids <- t(test_stage-predictions)

rmses <- apply(test_resids, 2, rmse)
maes <- apply(test_resids, 2, mae)
rmses_horizons <- apply(test_resids, 1, rmse)[c(1,24,28,72,96)]
mean(rmses)
sd(rmses)
mean(maes)
min(test_resids)
max(test_resids)
rmses_horizons

par(mfrow=c(3,4))
for(i in 1:12){
  pred_window <- (i+forecast_horizon):(i+forecast_horizon+forecast_horizon-1)
  plot(stage_test[pred_window]~dt_test[pred_window], type='l',lwd=3,col='lightgrey', main=paste0('Iter ',i))
  lines(predictions[i,]~dt_test[pred_window], col='blue')
}


highTideBin <- matrix(NA, nrow=96, ncol=2352)
for(t in 1:dim(test_resids)[2]){
  highTideBin[,t] <- (stage_test[t:(t+95)]>0)
}
highTideBin <- highTideBin[1:24,]
test_resids <- test_resids[1:24,]
ht_rmse_sum <- c()
for(j in 1:dim(test_resids)[2]){ht_rmse_sum <- c(ht_rmse_sum, sqrt(mean((test_resids[,j][highTideBin[,j]])^2)))}
print(mean(ht_rmse_sum)*100)
lt_rmse_sum <- c()
for(j in 1:dim(test_resids)[2]){lt_rmse_sum <- c(lt_rmse_sum, sqrt(mean((test_resids[,j][!highTideBin[,j]])^2)))}
print(mean(lt_rmse_sum)*100)

n_forecast=2205; horizon=(n_forecast):(n_forecast+95); i=1#2205 #420
rng=range(c(predictions[n_forecast,],stage_test[horizon]))

pdf(paste0('C:/Users/Jono/Dropbox/BackUp/E/Masters Advanced Analytics/Dissertation/msc_coastal_flooding/Images/96-hour Forecasts/','LSTM24-',n_forecast,'.pdf'), width = 9, height = 4)
par(mfrow=c(1,1), mar=c(4,4,2,0.5)); 
plot(stage_test[horizon]~dt_test[horizon], col='darkgrey', pch=19, cex=1.2,
     main=paste0('LSTM24',' 96-hr Forecast starting ',dt_test[horizon][1]), 
     ylab='Stage (m)',xlab='Datetime',
     ylim=rng)
lines(predictions[n_forecast,]~dt_test[horizon], col='blue', lwd=2)
dev.off()
