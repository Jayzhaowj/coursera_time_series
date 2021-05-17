######################
#####  Week 3  #######
######################

#################################################
##### Univariate DLM: Known, constant variance
#################################################

##### First order polynomial #####
### forward update equations ###
forward_filter <- function(data, parameters, initial_states){
  
  ## retrieve dataset
  yt <- data$yt
  num_timesteps <- length(yt)
  
  ## retrieve a set of quadruples 
  # FF, GG, VV, WW are scalar
  FF <- parameters$FF  
  GG <- parameters$GG
  VV <- parameters$VV
  WW <- parameters$WW
  
  ## retrieve initial states
  m0 <- initial_states$m0
  C0 <- initial_states$C0

  ## create placeholder for results
  d <- dim(GG)[1]
  at <- matrix(NA, nrow=num_timesteps, ncol=d)
  Rt <- array(NA, dim=c(d, d, num_timesteps))
  ft <- numeric(num_timesteps)
  Qt <- numeric(num_timesteps)
  mt <- matrix(NA, nrow=num_timesteps, ncol=d)
  Ct <- array(NA, dim=c(d, d, num_timesteps))
  et <- numeric(num_timesteps)
  
  
  for(i in 1:num_timesteps){
    # moments of priors at t
    if(i == 1){
      at[i, ] <- GG %*% t(m0)
      Rt[, , i] <- GG %*% C0 %*% t(GG) + WW
    }else{
      at[i, ] <- GG %*% t(mt[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , i-1] %*% t(GG) + WW
    }
    
    # moments of one-step forecast:
    ft[i] <- t(FF) %*% t(at[i, ]) 
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
    
    # moments of posterior at t:
    At <- Rt[, , i] %*% FF / Qt[i]
    et[i] <- yt[i] - ft[i]
    mt[i, ] <- at[i, ] + t(At) * et[i]
    Ct[, , i] <- Rt[, , i] - Qt[i] * At %*% t(At)
  }
  cat("Forward filtering is completed!") # indicator of completion
  return(list(mt = mt, Ct = Ct, at = at, Rt = Rt, ft = ft, Qt = Qt))
}

### smoothing function ###
backward_smoothing <- function(data, parameters, posterior_states){
  ## retrieve data 
  yt <- data$yt
  num_timesteps <- length(yt) 
    
  ## retrieve parameters
  FF <- parameters$FF
  GG <- parameters$GG
  
  ## retrieve parameters
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  at <- posterior_states$at
  Rt <- posterior_states$Rt
  
  ## create placeholder for posterior moments 
  mnt <- matrix(NA, nrow = dim(mt)[1], ncol = dim(mt)[2])
  Cnt <- array(NA, dim = dim(Ct))
  fnt <- numeric(num_timesteps)
  Qnt <- numeric(num_timesteps)
  for(i in num_timesteps:1){
    if(i == num_timesteps){
      mnt[i, ] <- mt[i, ]
      Cnt[, , i] <- Ct[, , i]
    }else{
      inv_Rtp1 <- chol2inv(chol(Rt[, , i+1]))
      Bt <- Ct[, , i] %*% GG %*% inv_Rtp1
      mnt[i, ] <- mt[i, ] + Bt %*% (mnt[i+1, ] - at[i+1, ])
      Cnt[, , i] <- Ct[, , i] + Bt %*% (Cnt[, , i + 1] - Rt[, , i+1]) %*% t(Bt)
    }
    fnt[i] <- t(FF) %*% t(mnt[i, , drop=FALSE])
    Qnt[i] <- t(FF) %*% t(Cnt[, , i]) %*% FF
  }
  cat("Backward smoothing is completed!")
  return(list(mnt = mnt, Cnt = Cnt, fnt=fnt, Qnt=Qnt))
}


## Forecast Distribution for k step
forecast_function <- function(posterior_states, k, parameters){

  ## retrieve parameters
  FF <- parameters$FF
  GG <- parameters$GG
  WW <- parameters$WW
  VV <- parameters$VV
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  
  ## set up parameters
  num_timesteps <- dim(mt)[1] # time points
  d <- dim(mt)[2]             # dimension of state parameters
  
  ## placeholder for results
  at <- matrix(NA, nrow = k, ncol = d)
  Rt <- array(NA, dim=c(d, d, k))
  ft <- numeric(k)
  Qt <- numeric(k)
  
  
  for(i in 1:k){
    ## moments of state distribution
    if(i == 1){
      at[i, ] <- GG %*% t(mt[num_timesteps, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , num_timesteps] %*% t(GG) + WW
    }else{
      at[i, ] <- GG %*% t(at[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Rt[, , i-1] %*% t(GG) + WW
    }
    
    ## moments of forecast distribution
    ft[i] <- t(FF) %*% t(at[i, , drop=FALSE])
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
  }
  cat("Forecasting is completed!") # indicator of completion
  return(list(at=at, Rt=Rt, ft=ft, Qt=Qt))
}

## create list for parameters
set_up_parameters <- function(GG, FF, VV, WW){
  return(list(GG=GG, FF=FF, VV=VV, WW=WW))
}


## create list for initial states
set_up_initial_states <- function(m0, C0){
  return(list(m0=m0, C0=C0))
}

## obtain 95% credible interval
get_credible_interval <- function(ft, Qt, quantile = c(0.025, 0.975)){
  z_quantile <- qnorm(quantile)
  bound <- matrix(0, nrow=length(ft), ncol=2)
  bound[, 1] <- ft + z_quantile[1]*sqrt(as.numeric(Qt)) # lower bound
  bound[, 2] <- ft + z_quantile[2]*sqrt(as.numeric(Qt)) # upper bound
  return(bound)
}

myplot <- function(results_filtered, results_smoothed, results_forecast, ci=TRUE){
  par(mfrow=c(1,1), cex.main=1.5, cex.lab=1.5)
  
  index <- seq(1871, 1970, length.out = length(Nile))
  index_train <- index[1:length(data$yt)]
  index_valid <- index[(length(data$yt)+1):length(Nile)]
  
  plot(index, Nile, pch=1, ylab = "level", main = "Nile River Level")
  lines(index_train, results_filtered$mt, type='l', col='red')
  lines(index_valid, results_forecast$ft, type='l', col='green')
  lines(index_train, results_smoothed$mnt, type='l', col='blue')
  if(ci){
    lines(index_train, ci_filtered[, 1], type='l', col='red', lty=2)
    lines(index_train, ci_filtered[, 2], type='l', col='red', lty=2)
    lines(index_valid, ci_forecast[, 1], type='l', col='green', lty=2)
    lines(index_valid, ci_forecast[, 2], type='l', col='green', lty=2)
    lines(index_train, ci_smoothed[, 1], type='l', col='blue', lty=2)
    lines(index_train, ci_smoothed[, 2], type='l', col='blue', lty=2)
  }
  legend('topright', legend=c("data", "filtered", "smoothed", "one-step-ahead"),
         col = c("black", "red", "blue", "green"), lty=c(1, 1, 1, 1))
}

##########################################################
#### first order polynomial with all parameters known ####
##########################################################
library(dlm)
num_total <- length(Nile)
num_timesteps <- 95
train_data <- Nile[1:num_timesteps]
valid_data <- Nile[(num_timesteps+1):num_total]

data <- list(yt = train_data)

## set up parameters
GG <- as.matrix(1)
FF <- as.matrix(1)
VV <- as.matrix(15100)
WW <- as.matrix(755)
m0 <- as.matrix(0.0)
C0 <- as.matrix(1e7)

## wrap up all parameters and initial values
initial_states <- set_up_initial_states(m0, C0)
parameters <- set_up_parameters(GG, FF, VV, WW)

## filtering
results_filtered <- forward_filter(data, parameters, initial_states)
ci_filtered <- get_credible_interval(results_filtered$mt, results_filtered$Ct)

## one step-ahead forecasting
results_forecast <- forecast_function(results_filtered, length(valid_data), parameters)
ci_forecast <- get_credible_interval(results_forecast$ft, results_forecast$Qt)

## smoothing
results_smoothed <- backward_smoothing(data, parameters, results_filtered)
ci_smoothed <- get_credible_interval(results_smoothed$mnt, results_smoothed$Cnt)


## plot results
myplot(results_filtered = results_filtered, 
       results_smoothed = results_smoothed, 
       results_forecast = results_forecast, ci=TRUE)


###############################################################
#### Univariate DLM: unknown, constant variance V = 1/phi #####
###############################################################
##### more generic version
### forward update equations ###
forward_filter <- function(data, parameters, initial_states, delta){
  
  ## retrieve datasets
  yt <- data$yt
  num_timesteps<- length(yt)
  
  ## retrieve parameters
  Ft <- parameters$Ft
  Gt <- parameters$Gt
  if(missing(delta)){
    Wt_star <- parameters$Wt_star
  }
  
  ## retrieve initial state
  m0 <- initial_states$m0
  C0_star <- initial_states$C0_star
  n0 <- initial_states$n0
  S0 <- initial_states$S0
  C0 <- S0*C0_star
  
  ## create placeholder for results
  d <- dim(Gt)[1]
  at <- matrix(NA, nrow=num_timesteps, ncol=d)
  Rt <- array(NA, dim=c(d, d, num_timesteps))
  ft <- numeric(num_timesteps)
  Qt <- numeric(num_timesteps)
  mt <- matrix(NA, nrow=num_timesteps, ncol=d)
  Ct <- array(NA, dim=c(d, d, num_timesteps))
  et <- numeric(num_timesteps)
  nt <- numeric(num_timesteps)
  St <- numeric(num_timesteps)
  dt <- numeric(num_timesteps)
  
  # moments of priors at t
  for(i in 1:num_timesteps){
    if(i == 1){
      at[i, ] <- Gt[, , i] %*% t(m0)
      Pt <- Gt[, , i] %*% C0 %*% t(Gt[, , i])
      if(missing(delta)){
        Wt <- Wt_star[, , i]*S0
        Rt[, , i] <- Pt + Wt
      }else{
        Rt[, , i] <- Pt/delta
      }
      
    }else{
      at[i, ] <- Gt[, , i] %*% t(mt[i-1, , drop=FALSE])
      Pt <- Gt[, , i] %*% Ct[, , i-1] %*% t(Gt[, , i])
      if(missing(delta)){
        Wt <- Wt_star[, , i] * St[i-1]
        Rt[, , i] <- Pt + Wt
      }else{
        Rt[, , i] <- Pt/delta
      }
    }
    
    # moments of one-step forecast:
    ft[i] <- t(Ft[, , i]) %*% t(at[i, , drop=FALSE]) 
    Qt[i] <- t(Ft[, , i]) %*% Rt[, , i] %*% Ft[, , i] + ifelse(i==1, S0, St[i-1])
    et[i] <- yt[i] - ft[i]
    
    nt[i] <- ifelse(i==1, n0, nt[i-1]) + 1
    St[i] <- ifelse(i==1, S0, St[i-1])*(1 + 1/nt[i]*(et[i]^2/Qt[i]-1))
   
    # moments of posterior at t:
    At <- Rt[, , i] %*% Ft[, , i] / Qt[i]

    mt[i, ] <- at[i, ] + t(At) * et[i]
    Ct[, , i] <- St[i]/ifelse(i==1, S0, St[i-1])*(Rt[, , i] - Qt[i] * At %*% t(At))
  }
  cat("Forward filtering is completed!\n")
  return(list(mt = mt, Ct = Ct, 
              at = at, Rt = Rt, 
              ft = ft, Qt = Qt, 
              nt = nt, St = St))
}

### smoothing function ###
backward_smoothing <- function(data, parameters, posterior_states, delta){
  ## retrieve data 
  yt <- data$yt
  num_timesteps <- length(yt) 
  
  ## retrieve parameters
  if(missing(delta)){
    Ft <- parameters$Ft
    Gt <- parameters$Gt
  }

  ## retrieve parameters
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  Rt <- posterior_states$Rt
  nt <- posterior_states$nt
  St <- posterior_states$St
  at <- posterior_states$at
  
  ## create placeholder for posterior moments 
  mnt <- matrix(NA, nrow = dim(mt)[1], ncol = dim(mt)[2])
  Cnt <- array(NA, dim = dim(Ct))
  fnt <- numeric(num_timesteps)
  Qnt <- numeric(num_timesteps)
  
  for(i in num_timesteps:1){
    if(i == num_timesteps){
      mnt[i, ] <- mt[i, ]
      Cnt[, , i] <- Ct[, , i]
    }else{
      if(missing(delta)){
        inv_Rtp1 <- chol2inv(chol(Rt[, , i+1]))
        Bt <- Ct[, , i] %*% Gt[, , i+1] %*% inv_Rtp1
        mnt[i, ] <- mt[i, ] + Bt %*% (mnt[i+1, ] - at[i+1, ])
        Cnt[, , i] <- Ct[, , i] + Bt %*% (Cnt[, , i+1] - Rt[, , i+1]) %*% t(Bt)
      }else{
        inv_Gt <- chol2inv(chol(Gt[, , i+1]))
        mnt[i, ] <- (1-delta)*mt[i, ] + delta*inv_Gt %*% t(mnt[i+1, , drop=FALSE])
        Cnt[, , i] <- St[num_timesteps]/St[i]*((1-delta)*Ct[, , i] + delta^2*inv_Gt %*% Cnt[, , i + 1]  %*% t(inv_Gt))
      }
    }
    fnt[i] <- t(Ft[, , i]) %*% t(mnt[i, , drop=FALSE])
    Qnt[i] <- t(Ft[, , i]) %*% t(Cnt[, , i]) %*% Ft[, , i]
  }
  cat("Backward smoothing is completed!\n")
  return(list(mnt = mnt, Cnt = Cnt, fnt=fnt, Qnt=Qnt))
}

## Forecast Distribution for k step
forecast_function <- function(posterior_states, k, parameters, delta){
  
  ## retrieve parameters
  Ft <- parameters$Ft
  Gt <- parameters$Gt
  if(missing(delta)){
    Wt_star <- parameters$Wt_star
  }
  
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  St <- posterior_states$St
  at <- posterior_states$at
  
  ## set up parameters
  num_timesteps <- dim(mt)[1] # time points
  d <- dim(mt)[2]             # dimension of state parameters
  
  ## placeholder for results
  at <- matrix(NA, nrow = k, ncol = d)
  Rt <- array(NA, dim=c(d, d, k))
  ft <- numeric(k)
  Qt <- numeric(k)
  
  for(i in 1:k){
    ## moments of state distribution
    if(i == 1){
      at[i, ] <- Gt[, , num_timesteps+i] %*% t(mt[num_timesteps, , drop=FALSE])
      
      if(missing(delta)){
        Rt[, , i] <- Gt[, , num_timesteps+i] %*% Ct[, , num_timesteps] %*% t(Gt[, , num_timesteps+i]) + St[num_timesteps]*Wt_star[, , num_timesteps+i]
      }else{
        Rt[, , i] <- Gt[, , num_timesteps+i] %*% Ct[, , num_timesteps] %*% t(Gt[, , num_timesteps+i])/delta
      }
      
    }else{
      at[i, ] <- Gt[, , num_timesteps+i] %*% t(at[i-1, , drop=FALSE])
      if(missing(delta)){
        Rt[, , i] <- Gt[, , num_timesteps+i] %*% Rt[, , i-1] %*% t(Gt[, , num_timesteps+i]) + St[num_timesteps]*Wt_star[, , num_timesteps + i]
      }else{
        Rt[, , i] <- Gt[, , num_timesteps+i] %*% Rt[, , i-1] %*% t(Gt[, , num_timesteps+i])/delta
      }
    }
    
    ## moments of forecast distribution
    ft[i] <- t(Ft[, , num_timesteps+i]) %*% t(at[i, , drop=FALSE])
    Qt[i] <- t(Ft[, , num_timesteps+i]) %*% Rt[, , i] %*% Ft[, , num_timesteps+i] + St[num_timesteps]
  }
  cat("Forecasting is completed!\n") # indicator of completion
  return(list(at=at, Rt=Rt, ft=ft, Qt=Qt))
}

## create list for parameters
set_up_parameters <- function(Gt, Ft, Wt_star){
  if(!is.array(Gt)){
    Stop("Gt and Ft should be array")
  }
  if(missing(Wt_star)){
    return(list(Gt=Gt, Ft=Ft))
  }else{
    return(list(Gt=Gt, Ft=Ft, Wt_star=Wt_star))
  }
}


## create list for initial states
set_up_initial_states <- function(m0, C0_star, n0, S0){
  return(list(m0=m0, C0_star=C0_star, n0=n0, S0=S0))
}

## obtain 95% credible interval
get_credible_interval <- function(ft, Qt, nt, quantile = c(0.025, 0.975)){
  bound <- matrix(0, nrow=length(ft), ncol=2)
  # lower bound of 95% credible interval
  z_quantile <- qt(quantile[1], df = nt)
  bound[, 1] <- ft + z_quantile*sqrt(as.numeric(Qt)) 
  
  # upper bound of 95% credible interval
  z_quantile <- qt(quantile[2], df = nt)
  bound[, 2] <- ft + z_quantile*sqrt(as.numeric(Qt)) 
  return(bound)
}

## separate the dataset into train data and valid data
num_total <- length(Nile)
num_timesteps <- 95
train_data <- Nile[1:num_timesteps]
valid_data <- Nile[(num_timesteps+1):num_total]
data <- list(yt = train_data)


## set up parameters
Gt <- array(1, dim = c(1, 1, num_total))
Ft <- array(1, dim = c(1, 1, num_total))
Wt_star <- array(100, dim = c(1, 1, num_total))
m0 <- as.matrix(0.0)
C0_star <- as.matrix(1e4)
n0 <- 10
S0 <- 1/10

## wrap up all parameters and initial values
initial_states <- set_up_initial_states(m0, C0_star, n0, S0)
parameters <- set_up_parameters(Gt, Ft, Wt_star)

## filtering
results_filtered <- forward_filter(data, parameters, initial_states)
ci_filtered <- get_credible_interval(results_filtered$mt, results_filtered$Ct, results_filtered$nt)

## smoothing
results_smoothed <- backward_smoothing(data, parameters, results_filtered)
ci_smoothed <- get_credible_interval(results_smoothed$mnt, results_smoothed$Cnt, 
                                     results_filtered$nt[num_timesteps])

## one-step ahead forecasting
results_forecast <- forecast_function(results_filtered, length(valid_data), parameters)
ci_forecast <- get_credible_interval(results_forecast$ft, results_forecast$Qt, 
                                     results_filtered$nt[num_timesteps])


## plot results
myplot(results_filtered = results_filtered,
       results_smoothed = results_smoothed,
       results_forecast = results_forecast)

##################################################
##### using discount factor ##########
##################################################

## compute log likelihood of one step ahead forecast function
log_likelihood_one_step_ahead <- function(yt, ft, Qt, nt){
  ## yt is time series, dimension: num_timesteps*1
  ## ft is expectation of one-step-ahead forecast function, dimension: num_timesteps*1
  ## Qt is variance of one-step-ahead forecast function, dimension: num_timesteps*1
  ## nt is degree freedom of t distribution, including initial states, dimension: (num_timesteps+1)*1
  num_timesteps <- length(yt)
  ntm1 <- nt[1:num_timesteps] # n_{t-1}
  zt <- (yt-ft)/sqrt(Qt) # standardization
  return(sum(dt(zt, df=ntm1, log=TRUE)))
}

## Maximize log density of one-step-ahead forecast function to select discount factor
adaptive_dlm <- function(data, parameters, initial_states, df_range){
  ll <- NA
  ll_tmp <- numeric(length(df_range))
  df_opt <- NA
  j <- 0
  ## find the optimal discount factor
  for(i in df_range){
    j <- j + 1
    results_tmp <- forward_filter(data, parameters, initial_states, i)
    ll_tmp[j] <- log_likelihood_one_step_ahead(data$yt, results_tmp$ft, results_tmp$Qt, 
                                               c(initial_states$n0, results_tmp$nt))
    if(j == 1){
      ll <- ll_tmp[j]
      results_filtered <- results_tmp
      df_opt <- i
    }else if(ll_tmp[j] > ll){
      ll <- ll_tmp[j]
      results_filtered <- results_tmp
      df_opt <- i
    }
  }
  results_smoothed <- backward_smoothing(data, parameters, results_filtered, delta = df_opt)
  results_forecast <- forecast_function(results_filtered, length(valid_data), parameters, df_opt)
  return(list(results_filtered=results_filtered, 
              results_smoothed=results_smoothed, 
              results_forecast=results_forecast, 
              df_opt = df_opt, ll = ll_tmp))
}


## separate the dataset into train data and valid data
num_total <- length(Nile)
num_timesteps <- 95
train_data <- Nile[1:num_timesteps]
valid_data <- Nile[(num_timesteps+1):num_total]
data <- list(yt = train_data)


## set up parameters
Gt <- array(1, dim = c(1, 1, num_total))
Ft <- array(1, dim = c(1, 1, num_total))
m0 <- as.matrix(0.0)
C0_star <- as.matrix(1e4)
n0 <- 10
S0 <- 1/10

## wrap up all parameters and initial values
initial_states <- set_up_initial_states(m0, C0_star, n0, S0)
parameters <- set_up_parameters(Gt, Ft)

## set up range of discount factor
df_range <- seq(0.8, 1, by = .01)

## fit discount DLM
results <- adaptive_dlm(data, parameters, initial_states, df_range)

## print selected discount factor
print(paste("The selected discount factor:", results$df_opt))

## retrieve filtered results
results_filtered <- results$results_filtered
ci_filtered <- get_credible_interval(results_filtered$mt, results_filtered$Ct, results_filtered$nt)

## retrieve smoothed results
results_smoothed <- results$results_smoothed
ci_smoothed <- get_credible_interval(results_smoothed$mnt, results_smoothed$Cnt, 
                                     results_filtered$nt[num_timesteps])

## retrieve one-step-ahead forecast results
results_forecast <- results$results_forecast
ci_forecast <- get_credible_interval(results_forecast$ft, results_forecast$Qt, results_filtered$nt[num_timesteps])

## plot results
myplot(results_filtered = results_filtered,
       results_smoothed = results_smoothed,
       results_forecast = results_forecast,
       ci = TRUE)
