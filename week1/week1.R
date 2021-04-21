###############################
#####       Week 1        #####
###############################


###############################
#####    stationarity     #####
###############################

#### simulation example: white noise 

## generate stationary time series
t <- 0:100
y_stationary <- rnorm(length(t), mean=0, sd=1) # the stationary time series (ts)

par(mfrow = c(1, 2), cex.lab = 1.3, cex.main = 1.3)
## plot the stationary signal and ACF
plot(t, y_stationary, type = 'l', col='red',
     xlab = 'time (t)', ylab = "Y(t)", 
     main = "Stationary signal")
acf(y_stationary, lag.max = length(y_stationary),
    xlab = "lag", ylab = "ACF", main = " ")


###########################################
##### differencing and moving average #####
###########################################

data(co2) # load co2 dataset
plot.ts(co2)

## doing moving average to remove effects of periodic component
mav_even <- function(x,d=12){filter(x,c(0.5/d,rep(1/d,(d-1)),0.5/d), sides=2)}
co2_mav <- mav_even(co2, d=12)

## plot smoothed time series
plot.ts(co2_mav, main = "Smoothed Series")

## plot sample ACF
acf(co2_mav[!is.na(co2_mav)], lag.max = 100, type = 'correlation', 
    main = "sample ACF of smoothed series")

## plot sample PACF
pacf(co2_mav[!is.na(co2_mav)], lag.max = 100, main = "sample PACF of smoothed series")

## doing 1st order differencing to remove effects of trend component
co2_mav_dtd <- diff(co2_mav, lag = 6, differences = 1)

## plot dataset after differencing
plot.ts(co2_mav_dtd, main = "Smoothed and Detrended Series")

## plot sample ACF
acf(co2_mav_dtd[!is.na(co2_mav_dtd)], lag.max = 50, type = "correlation", 
    main = "sample ACF of detrended and non-seasonal dataset")
## plot sample PACF
pacf(co2_mav_dtd[!is.na(co2_mav_dtd)], lag.max = 50, 
    main = "sample PACF of detrended and non-seasonal dataset")



##############################
#### only doing the difference ####
##############################

co2_dtd <- diff(diff(co2, lag=12, differences=1))

## plot sample ACF
acf(co2_dtd[!is.na(co2_dtd)], lag.max = 50, type = "correlation", 
    main = "sample ACF of detrended and non-seasonal dataset")
## plot sample PACF
pacf(co2_dtd[!is.na(co2_dtd)], lag.max = 50, 
     main = "sample PACF of detrended and non-seasonal dataset")




###############################
#####    sample ar(1)     #####
###############################
#### Example 1:
#### sample ar(1) with ar coefficients phi = 0.9 
#### and innovation standard deviation sigma = 0.5

set.seed(2021)
phi <- 0.9 # ar coefficient
sd <- 0.5 # innovation standard deviation
n <- 1000 # number of time points
lag.max <- 50 # max lag
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd) # generate stationary AR(1) process
par(mfrow = c(2, 2), cex.lab = 1.5)
ts.plot(yt) # plot 

## draw true ACF of AR(1)
cov_0 <- sd^2/(1-phi^2) # compute auto-covariance at h=0
cov_h <- phi^(0:lag.max)*cov_0 # compute auto-covariance at h != 0
plot(0:lag.max, cov_h/cov_0, pch = 1, type = 'h', col = 'red',
     ylab = "true ACF", xlab = "lag")


## draw sample autocorrelation function
acf(yt, lag.max = lag.max, 
    type = "correlation", ylab = "sample ACF", 
    lty = 1, ylim = c(0, 1), main = " ")

## draw sample partial autocorrelation function
pacf(yt, lag.ma = lag.max, main = "sample PACF")

#### Example 2:
#### sample ar(1) with ar coefficients phi = -0.9 
#### and innovation standard deviation sigma = 0.5
set.seed(2021)
phi <- -0.9 # ar coefficient
sd <- 0.5 # innovation standard deviation
n <- 1000 # number of time points
lag.max <- 50 # max lag
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd) # generate stationary AR(1) process
par(mfrow = c(2, 2), cex.lab = 1.5)
ts.plot(yt) # plot 

## draw true ACF of AR(1)
cov_0 <- sd^2/(1-phi^2) # compute auto-covariance at h=0
cov_h <- phi^(0:lag.max)*cov_0 # compute auto-covariance at h != 0
plot(0:lag.max, cov_h/cov_0, pch = 1, type = 'h', col = 'red',
     ylab = "true ACF", ylim = c(-1, 1), xlab = "lag")


## draw sample autocorrelation function
acf(yt, lag.max = lag.max, 
    type = "correlation", ylab = "sample ACF", 
    lty = 1, ylim = c(-1, 1), main = " ")

## draw sample partial autocorrelation function
pacf(yt, lag.ma = lag.max, main = "sample PACF")

#######################################################
#####           general AR(p) process           #######
#######################################################


########################################################
##### estimates of reciprocal characteristic roots #####
########################################################
## given AR coefficients
phi <- c(0.27, 0.07, -0.13, -0.15, -0.11, -0.15, -0.23, -0.14)
roots <- polyroot(c(1, -phi)) # compute reciprocal characteristic roots
modulus <- Mod(1/roots) # compute modulus
argument <- 2*pi/Arg(1/roots) # compute frequency
# print results modulus and frequency by decreasing order
print(cbind(modulus, argument)[order(modulus, decreasing=TRUE), ]) 

## given modulus and frequency, compute AR coefficients
modulus <- c(0.4898979, 0.4898979) # modulus
argument <- c(-6.891435, 6.891435) # frequency
phi <- numeric(2)
phi[1] <- Re(sum(modulus*exp(2*pi/argument*1i)))
phi[2] <- -Re(prod(modulus*exp(2*pi/argument*1i)))

## generate AR(2)
n <- 1000 # number of time points
sd <- 0.5 # innovation standard deviation
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd)

par(mfrow = c(1, 3), cex.lab = 1.5)
ts.plot(yt) # plot 


## draw sample autocorrelation function
acf(yt, lag.max = lag.max, 
    type = "correlation", ylab = "sample ACF", 
    lty = 1, ylim = c(-1, 1), main = " ")

## draw sample partial autocorrelation function
pacf(yt, lag.ma = lag.max, main = "sample PACF")

