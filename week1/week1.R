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
co2_mav_dtd <- diff(co2_mav, lag = 1, differences = 1)

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

####################################################
#####             MLE for AR(1)               ######
####################################################
set.seed(2021)
phi <- 0.9 # ar coefficient
sd <- 1 # innovation standard deviation
n <- 1000 # number of time points
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd) # generate stationary AR(1) process

## Case 1: Conditional likelihood
y <- as.matrix(yt[2:n]) # response
X <- as.matrix(yt[1:(n-1)]) # design matrix
phi_MLE <- as.numeric((t(X)%*%y)/sum(X^2)) # MLE for phi
s2 <- sum((y - phi_MLE*X)^2)/(length(y) - 1) # MLE for v 

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "MLE of conditional likelihood for v: ", s2, "\n")



## Case 2: Full likelihood
# log likelihood function
log_p <- function(phi, yt){
    0.5*(log(1-phi^2) - sum((yt[2:n] - phi*yt[1:(n-1)])^2) - yt[1]^2*(1-phi^2))
}

# Use built-in optimization method to obtain maximum likelihood estimates
result <- optimize(log_p, c(-1, 1), tol = 0.0001, maximum = TRUE, yt = yt)
cat("\n MLE of full likelihood for phi: ", result$maximum)


#######################################################
######          direct sampling                 #######
#######################################################

# direct sample from conditional likelihood and reference prior
n_sample <- 1000 # sample size

## step 1: sample posterior distribution of v from inverse gamma distribution
nu_sample <- 1/rgamma(n_sample, (n-2)/2, sum((yt[2:n] - phi_MLE*yt[1:(n-1)])^2)/2)

## step 2: sample posterior distribution of phi from normal distribution
phi_sample <- rnorm(n_sample, mean = phi_MLE, sd = sqrt(nu_sample)/sum((yt[1:(n-1)])^2))

## plot histogram of posterior samples of phi and nu
par(mfrow = c(1, 2), cex.lab = 1.3)
hist(phi_sample, xlab = bquote(phi), 
     main = bquote("Histogram of "~phi), xlim = c(0.895, .915))
abline(v = phi, col = 'red')
hist(nu_sample, xlab = bquote(nu), main = bquote("Histogram of "~nu))
abline(v = sd, col = 'red')

#######################################################
########        MCMC sampling                 #########
#######################################################

# compute Q star
Q_star <- function(phi, y){
    n <- length(y)
    return(y[1]^2*(1-phi^2) + sum((y[2:n] - phi*y[1:(n-1)])^2))
}

# transformation of phi
transform <- function(phi){
    return(log((1-phi)/(phi+1)))
}
# inverse transformation of phi
inv_transform <- function(eta){
    return((1-exp(eta))/(1+exp(eta)))
}
# compute likelihood 
log_likl <- function(phi, eta, v, y){
    return(log(1-phi^2) - 0.5*Q_star(phi, y)/v + eta - 2*log(1+exp(eta)))
}

# metropolis hasting
mh <- function(eta_cur, v, step_size, y){
    acc <- 0
    phi_cur <- inv_transform(eta_cur)
    eta_cand <- rnorm(1, mean = eta_cur, sd = step_size)
    phi_cand <- inv_transform(eta_cand)
    accept.prob <- log_likl(phi_cand, eta_cand, v, y) - log_likl(phi_cur, eta_cur, v, y)
    
    if(log(runif(1)) < accept.prob){
        eta_cur <- eta_cand
        phi_cur <- phi_cand
        acc <- 1
    }
    return(list(eta_cur = eta_cur,
                phi_cur = phi_cur,
                acc = acc))
        
}

## auto tuning the stepwise
auto_tuning <- function(acc, iter, tuning, target_ar = 0.5){
    index <- (iter - 49):(iter)
    accept_rate <- mean(acc[index])
    n <- iter / 50
    if(accept_rate > target_ar){
        new_tuning <- tuning * exp(sqrt(1/n))
    }else{
        new_tuning <- tuning / exp(sqrt(1/n))
    }
    return(new_tuning)
}

## run metropolis hasting
MCMC_sample <- function(n_sample, burnin, y){
    ## storage of samples and acceptance
    phi_sample <- numeric(n_sample)
    v_sample <- numeric(n_sample)
    acc <- numeric(burnin + n_sample)
    ## initialization and set up step size
    phi_cur <- 0
    eta_cur <- transform(phi_cur)
    step_size <- 1
    for(i in 1:(burnin+n_sample)){
        # sample variance v
        v_cur <- 1/rgamma(n = 1, shape = n/2, rate = Q_star(phi_cur, y)/2)
        
        # sample location eta
        tmp <- mh(eta_cur = eta_cur, v = v_cur, step_size = step_size, y)
        
        # retrieve results
        phi_cur <- tmp$phi_cur
        eta_cur <- tmp$eta_cur
        acc[i] <- tmp$acc
        
        # tune step size
        if(i %% 50 == 0){
            step <- auto_tuning(acc = acc, iter = i, 
                                tuning = step_size)
        }
        
        # store the sample after burnin period
        if(i > burnin){
            phi_sample[i-burnin] <- phi_cur
            v_sample[i-burnin] <- v_cur
        }
    }
    return(list(phi_sample = phi_sample, 
                v_sample=v_sample, 
                acc_rate = mean(acc[1:n_sample+burnin])))
}

# run MCMC sampling with sample size = 10000 and burnin = 5000
result <- MCMC_sample(n_sample = 10000, burnin = 5000, y = yt)

## print Acceptance rate
cat("\n Acceptance Rates: ", result$acc_rate, "\n")
## plot results of phi
par(mfrow = c(2, 2), cex.lab = 1.3)
plot(result$phi_sample, ylab = "values", 
     xlab = "iterations", type = 'l', main = bquote("Traceplot of "~phi))
hist(result$phi_sample, xlab = bquote(phi), main = bquote("Histogram of "~phi))
abline(v=phi, col="red")
## plot results of v
plot(result$v_sample, ylab = "values",
     xlab = "iterations", type = 'l', main = "Traceplot of v")
hist(result$v_sample, xlab = "v", main = "Histogram of v")
abline(v=sd, col="red")

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

###################################################
