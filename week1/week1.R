###############################
#####       Week 1        #####
###############################


###############################
#####    stationarity     #####
###############################
#### simulation example: white noise 
### generate stationary time series
t <- 0:100
y_stationary <- rnorm(length(t), mean=0, sd=1) # the stationary time series (ts)

par(mfrow = c(1, 2), cex.lab = 1.3, cex.main = 1.3)
# plot the stationary signal and ACF
plot(t, y_stationary, type = 'l', col='red',
     xlab = 'time (t)', ylab = "Y(t)", 
     main = "Stationary signal")
acf(y_stationary, lag.max = length(y_stationary),
    xlab = "lag", ylab = "ACF", main = " ")


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

