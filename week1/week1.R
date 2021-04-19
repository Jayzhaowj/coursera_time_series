###############################
#####       Week 1        #####
###############################



###############################
#####    sample ar(1)     #####
###############################
set.seed(2021)
phi <- 0.9 # ar coefficient
sd <- 0.5 # innovation standard deviation
n <- 1000 # number of time points
lag.max <- 50 # max lag
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd) # generate stationary AR(1) process
plot(yt) # plot 

acf(yt, lag.max = lag.max, type = "correlation", ylim = c(-1, 1))


###############################
#####    stationarity     #####
###############################

### generate stationary time series
t <- 0:100
y_stationary <- rnorm(length(t), mean=1, sd=1) # the stationary time series (ts)

par(mfrow = c(1, 2), cex.lab = 1.3, cex.main = 1.3)
# plot the stationary signal and ACF
plot(t, y_stationary, type = 'l', col='red',
     xlab = 'time (t)', ylab = "Y(t)", 
     main = "Stationary signal")
acf(y_stationary, lag.max = length(y_stationary),
    xlab = "lag", ylab = "ACF", main = " ")
