#############################
###### quiz 1
#############################
dir <- "C:/Users/johnn/Documents/Courses/coursera_time_series/week1/quiz/"


# ts 1
set.seed(1)
n_t <- 100
beta_0 <- .5
beta_1 <- .1
et <- rnorm(n_t)
yt <- ts(beta_0 + beta_1 * (1:n_t) + et)
png(filename = paste0(dir, "ts1.png"))
plot(yt, ylab = "value")
dev.off()

# ts 2
yt <- arima.sim(list(ar = 0.8, sd =1), n=n_t)
png(filename = paste0(dir, "ts2.png"))
plot(yt, ylab = "value")
dev.off()

# ts 3
yt <- numeric(n_t)
for(i in 1:n_t){
  yt[i] <- sum(rnorm(sum(1:i)))
}

png(filename = paste0(dir, "ts3.png"))
plot(ts(yt), ylab = "value")
dev.off()











###################################################
lag.max <- 10
phi <- 0.8*c(1, rep(0,9))
png(filename = paste0(dir, "pacf1.png"))
plot(1:lag.max, phi, pch= 1, type = 'h', col = 'black', ylab = "PACF", xlab = "lag")
dev.off()

sd <- 1
phi <- 0.9
## draw true ACF of AR(1)
cov_0 <- sd^2/(1-phi^2) # compute auto-covariance at h=0
cov_h <- phi^(1:lag.max)*cov_0 # compute auto-covariance at h != 0
png(filename = paste0(dir, "pacf2.png"))
plot(1:lag.max, cov_h/cov_0, pch = 1, type = 'h', col = 'black',
     ylab = "PACF", xlab = "lag")
dev.off()



sd <- 1
phi <- -0.9
## draw true ACF of AR(1)
cov_0 <- sd^2/(1-phi^2) # compute auto-covariance at h=0
cov_h <- phi^(1:lag.max)*cov_0 # compute auto-covariance at h != 0
png(filename = paste0(dir, "pacf3.png"))
plot(1:lag.max, cov_h/cov_0, pch = 1, type = 'h', col = 'black',
     ylab = "PACF", xlab = "lag")
dev.off()
