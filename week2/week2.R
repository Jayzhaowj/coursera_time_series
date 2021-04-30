###############################
#####       Week 2      #####
###############################


###############################
#####     R  examples     #####
###############################
set.seed(2021)

### Example 1
### Simulate from AR(2) with two real reciprocal roots (e.g., 0.95 and 0.5)
recip_roots <- c(0.95, 0.5) ## two different real reciprocal roots
phi <- c(sum(recip_roots), -prod(recip_roots)) ## compute ar coefficients
n <- 1000 ## set up number of time points
sd <- 1 ## set up standard deviation
y <- arima.sim(model = list(ar=phi), n = n, sd = sd) # generate ar(2)


par(mfrow = c(1, 2), cex.lab = 1.5, cex.main = 1.5)
### plot sample ACF
acf(y, lag.max = 50, type = "correlation", 
    main = "sample ACF", )


### plot sample PACF
pacf(y, lag.max = 50, main = "sample PACF")


### Example 2
### Simulate from AR(2) with a pair of complex roots (e.g., r=0.95 and lambda = 12)
set.seed(2021)

r <- c(0.95, 0.95)
lambda <- c(-15, 15)

# placeholder for phi
phi <- numeric(2) 
phi[1] <- Re(sum(r*exp(2*pi/lambda*1i))) # ar coefficients at lag 1
phi[2] <- -Re(prod(r*exp(2*pi/lambda*1i))) # ar coefficients at lag 2

n <- 1000 # number of time points
sd <- 1 # standard deviation
y <- arima.sim(model = list(ar=phi), n = n, sd = sd) # generate ar(2)

par(mfrow = c(1, 2), cex.lab = 1.5, cex.main = 1.5)
### plot sample ACF
acf(y, lag.max = 50, type = "correlation", 
    main = "sample ACF", )


### plot sample PACF
pacf(y, lag.max = 50, main = "sample PACF")


set.seed(2021)
### Simulate from AR(3) with one real root 
### and a pair of complex roots (e.g., r=0.95 and lambda = 12 and real root with
### 0.8 modulus)
r <- c(0.95, 0.95, 0.8) ## modulus
lambda <- c(-12, 12) ## lambda
recip_roots <- c(r[1:2]*exp(2*pi/lambda*1i), r[3]) ## reciprocal roots
phi <- numeric(3) # placeholder for phi
phi[1] <- Re(sum(recip_roots)) # ar coefficients at lag 1
phi[2] <- -Re(recip_roots[1]*recip_roots[2] + recip_roots[1]*recip_roots[3] + recip_roots[2]*recip_roots[3]) # ar coefficients at lag 2
phi[3] <- Re(prod(recip_roots))
n <- 1000 # number of time points
sd <- 1 # standard deviation
y <- arima.sim(model = list(ar=phi), n = n, sd = sd) # generate ar(2)

par(mfrow = c(1, 2), cex.lab = 1.5, cex.main = 1.5)
### plot sample ACF
acf(y, lag.max = 50, type = "correlation", 
    main = "sample ACF", )


### plot sample PACF
pacf(y, lag.max = 50, main = "sample PACF")



################################################################
####################################################
#####             MLE for AR(2)               ######
####################################################
set.seed(2021)
phi <- c(0.9, -0.2) # ar coefficient
sd <- 1 # innovation standard deviation
n <- 1000 # number of time points
yt <- arima.sim(n = n, model = list(ar = phi), sd = sd) # generate stationary AR(1) process

## Case 1: Conditional likelihood
p <- 2
y <- as.matrix(yt[(p+1):n]) # response
X <- cbind(yt[(p+1):n-1], yt[(p+1):n-2]) # design matrix
XtX <- t(X)%*%X
XtX_inv <- solve(XtX)
phi_MLE <- XtX_inv%*%t(X)%*%y # MLE for phi
s2 <- sum((y - X%*%phi_MLE)^2)/(length(y) - p) # MLE for v 

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "MLE of conditional likelihood for v: ", s2, "\n")



#######################################################
######          direct sampling                 #######
#######################################################

# direct sample from conditional likelihood and reference prior
n_sample <- 1000 # sample size

## step 1: sample posterior distribution of v from inverse gamma distribution
nu_sample <- 1/rgamma(n_sample, (n-1-p)/2, sum((y-X%*%phi_MLE)^2)/2)

## step 2: sample posterior distribution of phi from normal distribution
phi_sample <- matrix(NA, nrow = n_sample, ncol = p)
for(i in 1:n_sample){
    z <- rnorm(p)
    A <- t(chol(nu_sample[i]*XtX_inv))
    phi_sample[i, ] <- phi_MLE + A%*%as.matrix(z)
}

## plot histogram of posterior samples of phi and nu
par(mfrow = c(1, 3), cex.lab = 1.3)
for(i in 1:2){
    hist(phi_sample[, i], xlab = bquote(phi), 
         main = bquote("Histogram of "~phi[.(i)]))
    abline(v = phi[i], col = 'red')
}

hist(nu_sample, xlab = bquote(nu), main = bquote("Histogram of "~nu))
abline(v = sd, col = 'red')






###########################################
######   compute AIC and BIC for eeg data   #######
###########################################
dir <- "/Users/johnn/Documents/Courses/coursera_time_series/week2/"
setwd(paste0(dir, "datasets"))
data <- scan("eeg")
data <- data - mean(data) #subtract its mean
n_Total <- length(data) # number of total time points
pmax <- 25 # the maximum of model order
xsave<-t(matrix(data[rev(rep((1:pmax),n_Total-pmax)+rep((0:(n_Total-pmax-1)),rep(pmax,n_Total-pmax)))],pmax,n_Total-pmax));
y <- rev(data[(pmax+1):n_Total])
n_cond <- length(y) # number of total time points - the maximum of model order

## compute MLE
my_MLE <- function(y, xsave, p){
    n <- length(y)
    x <- xsave[,1:p]
    a <- solve(t(x) %*%x)
    #a <- (a + t(a))/2 # numerically stable
    b <- a%*%t(x)%*%y # mle for ar coefficients
    r <- y - x%*%b # residuals 
    nu <- n - p # degree freedom
    R <- sum(r*r) # SSE
    s <- R/nu #MSE
    return(list(b = b, s = s, R = R, nu = nu))
}


## function of AIC and BIC computation 
AIC_BIC <- function(y, xsave, p){
    ## number of time points
    n <- length(y)
    
    ## compute MLE
    tmp <- my_MLE(y, xsave, p)
    
    ## retrieve results
    R <- tmp$R
    
    ## compute maximum likelihood
    maxlikl <- n*log(R)
    
    ## compute AIC and BIC
    aic <- maxlikl + 2*(p)
    bic <- maxlikl + log(n)*(p)
    return(list(aic = aic, bic = bic))
}
# Compute AIC, BIC, and the marginal likelihood
aic <- numeric(pmax)
bic <- numeric(pmax)

for(p in 1:pmax){
    tmp <- AIC_BIC(y, xsave, p)
    aic[p] <- tmp$aic
    bic[p] <- tmp$bic
    print(c(p, aic[p], bic[p])) # print AIC and BIC by model order
}

## compute difference between the value and its minimum
aic <- aic-min(aic) 
bic <- bic-min(bic) 

## draw plot of AIC, BIC, and the marginal likelihood
par(mfrow = c(1, 1))
matplot(1:pmax,matrix(c(aic,bic),pmax,2),ylab='value',
        xlab='AR order p',pch="ab", col = 'black', main = "AIC and BIC")
text(which.min(aic), aic[which.min(aic)], "a", col = 'red')
text(which.min(bic), bic[which.min(bic)], "b", col = 'red')

########################################################
p <- which.min(bic) # We set up moder order
print(paste0("The chosen model order: ", p))

##
y <- as.matrix(rev(data[(p+1):n_Total]))
x <- t(matrix(data[rev(rep((1:p),n_Total-p)+ rep((0:(n_Total-p-1)),rep(p,n_Total-p)))],p,n_Total-p));
tmp_MLE <- my_MLE(y, x, p)
## posterior mean from AR reference analysis
phi <- tmp_MLE$b ## posterior mean of AR coefficients
s <- sqrt(tmp_MLE$s) ## posterior mean of standard deviation

########################################################
##### estimates of reciprocal characteristic roots #####
########################################################
roots <- polyroot(c(1, -phi)) # compute reciprocal characteristic roots
modulus <- Mod(1/roots) # compute modulus
argument <- 2*pi/Arg(1/roots) # compute frequency
# print results modulus and frequency by decreasing order
print(cbind(modulus, argument)[order(modulus, decreasing=TRUE), ]) 


#############################################
#######    draw spectral density    #########
#############################################

### using spec.ar to draw spectral density
spec.ar(as.numeric(y), order = p, method = "mle",
        main = "EEG")

### using arma.spec to draw spectral density
#install.packages("astsa")
require("astsa")
library("astsa")
par(mfrow=c(3,1))
arma.spec(log="no", main="White Noise")
arma.spec(ma=.5, log="no", main="Moving Average")
arma.spec(ar=c(1,-.9), log="no", main="Autoregression")

## plot spectral density of EEG data with estimated ar coefficients
## and estimated innvovation variance
par(mfrow = c(1, 1))
arma.spec(ar=phi, var.noise = s^2, main = 'EEG')


### Remark: The difference between function spec.ar and arma.spec
### the function spec.ar: you pass the data into spec.ar, then the function will 
### fit with built-in ar function and draw spectral density based on estimations

### the function arma.spec: the function draw sepctral density based on the coefficients you specified
##########################################################


###########################################
##### differencing and moving average #####
###########################################


data(co2) # load co2 dataset
data <- ts(co2, start=c(1959, 1), frequency = 12)
par(mfrow = c(1, 3), cex.lab = 1.5)
plot.ts(data) # plot time series


##############################
#### only doing the difference ####
##############################

co2_dtd <- diff(diff(data, lag=12, differences=1))
co2_dtd <- co2_dtd[!is.na(co2_dtd)]
## plot sample ACF
acf(co2_dtd, lag.max = 50, type = "correlation", 
    main = "sample ACF of detrended and non-seasonal dataset")
## plot sample PACF
pacf(co2_dtd, lag.max = 50, 
     main = "sample PACF of detrended and non-seasonal dataset")

#################################
#####    fit ARIMA model    #####
#################################
require(forecast)
library(forecast)

## Based on the ACF and PACF above, 
## The significant spike at lag 1 in the ACF suggests a non-seasonal MA(1) component, 
## and the significant spike at lag 12 in the ACF suggests a seasonal MA(1) component.


# We try ARIMA(0, 1, 1)(0, 1, 1)_12 
fit <- Arima(data, order = c(0, 1, 1), 
             seasonal = list(order=c(0, 1, 1), period = 12))

## check the residuals
checkresiduals(fit, lag=48)

