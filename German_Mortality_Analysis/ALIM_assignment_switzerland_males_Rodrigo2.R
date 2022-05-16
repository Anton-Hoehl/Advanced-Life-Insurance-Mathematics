# ALIM assignment - Switzerland_males

library(tidyverse)
library(demography)
library(forecast)
library(dplyr)

##----input data from HMD.org---------------------------------------------------------------------------

Switzerland_males_1876_2020 <-
  read_table(file = "./6_2_data/mltper_1x1.txt",
             col_names = T,
             skip = 1)

Switzerland_exposures_1876_2020 <-
  read_table(file = "./6_2_data/Exposures_1x1.txt",
             col_names = T,
             skip = 1)

##----data prep----------------------------------------------------------------------

Switzerland_males_1970_2019_young <-
  Switzerland_males_1876_2020 %>%
  mutate(
    expo = Switzerland_exposures_1876_2020$Male,
    Age = replace(Age, Age == "110+", 110),
    #Transforming 110+ into 110
    Age = as.numeric(Age)
  ) %>%
  mutate(dexpo = expo * mx) %>%
  filter(Year >= 1970 &
           Year <= 2019) %>%               #Subset from year 1980 until most recent year (2019) (40 years in total)
  filter(Age <= 79)                          #Subset from age 0 until 89 because the remaining ages will be closed with Kannisto

Switzerland_males_1970_2019_old_age <-
  Switzerland_males_1876_2020 %>%
  mutate(
    expo = Switzerland_exposures_1876_2020$Male,
    Age = replace(Age, Age == "110+", 110),
    #Transforming 110+ into 110
    Age = as.numeric(Age)
  ) %>%
  mutate(dexpo = expo * mx) %>%
  filter(Year >= 1970 &
           Year <= 2019) %>%               #Subset from year 1980 until most recent year (2019) (40 years in total)
  filter(Age > 79)                             #old ages will be closed with Kannisto


# Transform the 3 observations with "dxt = 0" into "dxt = 1"
Switzerland_males_1970_2019_young["dexpo"][Switzerland_males_1970_2019_young["dexpo"] == 0] <- 1

years_swiss = 1970:max(Switzerland_males_1970_2019_young$Year) # May be needed later
ages = 0:79 #May be needed later
abc = expand.grid(Year = years_swiss, Age = ages) # May be needed later

## Plot of the log of the central death rate - we observe an improvement in the central death rates throughout the years

p_swiss = ggplot(Switzerland_males_1970_2019_young,
                 aes(x = Age, y = log(mx), group = Year)) +
  geom_line(aes(colour = Year), size = 1, linetype = 1) +
  scale_colour_gradientn(colours = rainbow(10)) +
  scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
  theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") + ggtitle("Switzerland - mx evolution")
p_swiss

# Calibration of the Poisson Likelihood (optimization of the Poisson Likelihood
#with univariate Newton-Raphson steps)



# ---------------------------------//-------------------------------
## llmaxM2B (Borrowed function because I couldnt use from the fitModels script)

llmaxM2B=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b3,k2,g3 are given
  #   solve for b2
  b21=b2
  b20=b21-1
  thetat=k2*ev*exp(b1+b3*g3)
  s1=sum(dv*k2*wv)
  while(abs(b21-b20) > 0.1)
  {
    b20=b21
    f0=sum((exp(b20*k2)*thetat)*wv)-s1
    df0=sum((exp(b20*k2)*k2*thetat)*wv)
    b21=b20-f0/df0
  }
  b21
}

## --------------------------------------//-----------------------------

## llmaxM2D (Borrowed function because I couldnt use from the fitModels script)
# IMP: THIS FUNCTION CONTAINS THE NEWTON RAPHSON STEPS AS WE DISCUSSED ON THE LECTURE SHEETS

llmaxM2D=function(b1,b2,b3,k2,g3,dv,ev,wv=1){
  #   b1,b2,b3,g3 are given
  #   solve for k2
  k21=k2
  k20=k21-1
  thetat=b2*ev*exp(b1+b3*g3)
  s1=sum(dv*b2*wv)
  while(abs(k21-k20) > 0.1)
  {
    k20=k21
    f0=sum((exp(k20*b2)*thetat)*wv)-s1
    df0=sum((exp(k20*b2)*b2*thetat)*wv)
    k21=k20-f0/df0
  }
  k21
}

## --------------------------------------//-----------------------------

# 1st method : using fit701 function
fit701=function(xv,yv,etx,dtx,wa){
  # Model M1
  # Lee-Carter model
  # log m(t,x) = beta1(x) + beta2(x).kappa2(t) + Poisson error
  # Inputs:
  #   xv = vector of ages, length n
  #   yv = vector of years, length m
  #   etx = m x n matrix of exposures
  #   dtx = m x n matrix of deaths
  #   wa = m x n matrix of weights (0 or 1)
  xv<-as.vector(unlist(xv))
  yv<-as.vector(unlist(yv))
  etx<-as.matrix(etx)
  dtx<-as.matrix(dtx)
  wa<-as.matrix(wa)
  mtx=dtx/etx	  # matrix of death rates
  
  qtx=1-exp(-mtx) # matrix of mortality rates
  
  if(max(xv) > 89)
  {  
    cat("Upper age too high - suggest abort programme\n") 
  }
  
  n=length(xv)	# number of ages
  m=length(yv)	# number of years
  
  cy=(yv[1]-xv[n]):(yv[m]-xv[1])  # cohort approximate years of birth
  
  # initialise parameter vectors
  beta1v=(1:n)*0
  beta2v=(1:n)*0
  beta3v=(1:n)*0		# dummy vector, this will stay at 0
  kappa2v=(1:m)*0
  gamma3v=(1:(n+m-1))*0	# dummy vector, this will stay at 0
  ia=array((1:m),c(m,n))	# matrix of year indexes, i, for the data
  ja=t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya=ia-ja		 	# matrix of year of birth indexes for the data
  imj=(1-n):(m-1)		# the range of values taken by i-j
  lg=n+m-1		 	# number of different values taken by i-j
  ca=ya+yv[1]-xv[1]		# matrix of years of birth
  
  # Now set weights to zero for cohorts with fewer than 5 observations
  for(k in 1:lg)
  {
    nk=sum((ca == cy[k])*wa)
    if(nk < 5)
    {
      wa=wa*(1- (ca == cy[k]))
    }
  }
  
  ww=cy*0+1	 # this is a vector of 1's and 0's with
  # a 0 if the cohort is completely excluded
  for(k in 1:lg)
  {
    ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
  }
  
  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
  mx=mean(xv)
  for(j in 1:n)
  {
    beta1v[j]=sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
    beta2v[j]=1/n
  }
  kappa2v=(m:1)-(m+1)/2
  
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  while(abs(l1-l0) > 0.0001)
  {
    iteration=iteration+1
    
    l0=l1
    # Stage 1B optimise over the beta2(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta2v[j]=llmaxM2B(beta1v[j],beta2v[j],beta3v[j],
                         kappa2v,gamma3v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1,"-> ")
    
    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]	# exposure
      kappa2v[i]=llmaxM2D(beta1v,beta2v,beta3v,
                          kappa2v[i],gamma3v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    # Now apply the constraints (IDENTIFIABILITY CONSTRAINTS)
    fac21=mean(kappa2v)
    fac22=sum(beta2v)
    kappa2v=fac22*(kappa2v-fac21)    # ensures that the kappas sum to 0
    beta2v=beta2v/fac22		     # ensures that the beta2's sum to 1
    beta1v=beta1v+beta2v*fac22*fac21 # => beta1 needs to be adjusted to compensate
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1," ->")
    
    # Stage 1A optimise over the beta1(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      wv=1	    # can be set to a vector of weights
      # to e.g. exclude duff years
      wv=wa[,j]
      s1=sum(wv*dtx[,j])
      s2=sum(wv*etx[,j]*exp(beta2v[j]*kappa2v+beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
      beta1v[j]=log(s1)-log(s2)
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v+beta2v*kappa2v[i]+beta3v*gamma3v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa)
    cat(l1,"\n")
    
  }		 # end while loop
  
  # calculate number of parameters and deduct 4 for the number of constraints
  npar=length(beta1v)+length(beta2v)+length(kappa2v)-2
  
  # Calculate the BIC
  BIC=l1-0.5*log(sum(wa))*npar
  
  list(beta1=beta1v,beta2=beta2v,beta3=beta3v,
       kappa2=kappa2v,gamma3=gamma3v,x=xv,y=yv,cy=cy,
       wa=wa,epsilon=epsilon,mhat=mhat,ll=l1,BIC=BIC,npar=npar,mtxLastYear = mtx[m,])		 
}

## --------------------------------------//-----------------------------
# Fit701 function arguments
#   xv = vector of ages, length n
#   yv = vector of years, length m
#   etx = m x n matrix of exposures
#   dtx = m x n matrix of deaths
#   wa = m x n matrix of weights (0 or 1)

ages_swiss = as.matrix(ages, nrow = 1)
years_swiss = as.matrix(years_swiss, nrow = 1)
ext_swiss = matrix(Switzerland_males_1970_2019_young$expo, nrow = 50, byrow = TRUE)
dxt_swiss = matrix(Switzerland_males_1970_2019_young$dexpo, nrow = 50, byrow = TRUE)

# Estimates for Bx(1), Bx(2), Kt(2)

LCfit701 <-  fit701(ages_swiss,years_swiss,ext_swiss,dxt_swiss,matrix(1, length(years_swiss), length(ages_swiss)))

names(LCfit701)

LCfit701$beta1 # gives the beta1 for each age (90 ages)
LCfit701$beta2 # gives the beta2 for each age (90 ages)
LCfit701$kappa2 # gives the kt for each year (50 years)

## We now check if the parameters beta2 and kappa2 satisfy the identifiability constraints:

sum(LCfit701$beta2) # = 1
sum(LCfit701$kappa2) # = 0 

##----------------------------------//-----------------------------------------

#### 2.2. Using self-written function ####
## Function ##
LCNRopt <- function(dxt, ext, eps = 1e-4, maxiter = 1e4) {
  mxt  = dxt / ext
  m    = ncol(ext)
  LL <- function(dxt, ext, Beta, Kappa) {
    Mat = matrix(NA, nrow(dxt), ncol(dxt))
    for(i in seq_len(nrow(dxt)))
      for(j in seq_len(ncol(dxt)))
        Mat[i, j] = dxt[i, j] * (Beta[i, 1] * Beta[i, 2] * Kappa[j]) - 
          ext[i, j] * exp(Beta[i, 1] * Beta[i, 2] * Kappa[j])
    sum(apply(Mat, 1, sum))
  }
  Beta  = cbind(apply(mxt, 1, function(x) sum(log(x))) / ncol(mxt), rep(1 / nrow(dxt), nrow(mxt)))
  Kappa = (m : 1) - (m + 1) / 2
  Conv  = F
  iter  = 0
  LogL  = NULL
  LogL[iter + 1] = LL(dxt, ext, Beta, Kappa)
  
  while(!Conv) {
    if((iter %% 1000) == 0)
      cat("\n\nIteration number", iter, "\n\n")
    for(i in seq_len(nrow(dxt))) {
      B0i = Beta[i, 1]
      B2i = Beta[i, 2]
      dxti = dxt[i, ]
      exti = ext[i, ]
      B1i = B0i - (sum(dxti - exti * exp(B0i + B2i * Kappa))) /
        - (sum(exti * exp(B0i + B2i * Kappa)))
      
      while(abs(B0i - B1i) > 0.01) {
        B0i = B1i
        B1i = B0i - (sum(dxti - exti * exp(B0i + B2i * Kappa))) /
          - (sum(exti * exp(B0i + B2i * Kappa)))
      }
      Beta[i, 1] = B1i
    }
    
    for(i in seq_len(ncol(dxt))) {
      B1 = Beta[, 1]
      B2 = Beta[, 2]
      dxti = dxt[, i]
      exti = ext[, i]
      K0i = Kappa[i] 
      K1i = K0i - (sum( (dxti -  exti * exp(B1 + B2 * K0i)) * B2 ))  / 
        - (sum(exti * exp(B1 + B2 * K0i) * B2^2))
      
      while(abs(K1i - K0i) > 0.01) {
        K0i = K1i
        K1i = K0i - (sum( (dxti -  exti * exp(B1 + B2 * K0i)) * B2 ))  /    ## Here we clearly see the Newton-Raphson step
          - (sum(exti * exp(B1 + B2 * K0i) * B2^2))                         ## for the Kt
      }
      Kappa[i] = K1i
    }
    
    SumB2     = sum(Beta[, 2])
    AvgKappa  = mean(Kappa)
    Kappa     = SumB2 * (Kappa - AvgKappa)
    Beta[, 2] = Beta[, 2] / SumB2
    Beta[, 1] = Beta[, 1] + Beta[, 2] * SumB2 * AvgKappa
    
    for(i in seq_len(nrow(dxt))) {
      B0i = Beta[i, 1]
      B2i = Beta[i, 2]
      dxti = dxt[i, ]
      exti = ext[i, ]
      
      B2.1i = B2i - (sum( (dxti - exti * exp(B0i + B2i * Kappa)) * Kappa )) /
        - (sum( exti * (exp(B0i + B2i * Kappa)) * Kappa^2) )
      
      while(abs(B2i - B2.1i) > 0.01) {
        B2i = B2.1i
        B2.1i = B2i - (sum( (dxti - exti * exp(B0i + B2i * Kappa)) * Kappa )) /
          - (sum( exti * (exp(B0i + B2i * Kappa)) * Kappa^2) )
      }
      
      Beta[i, 2] = B2.1i
    }
    
    iter = iter + 1 
    LogL[iter + 1] = LL(dxt, ext, Beta, Kappa)
    if(abs(LogL[iter + 1] - LogL[iter]) < eps)
      break
    if(iter > maxiter)
      break
  }
  if(iter > maxiter)
    warning("Maximum number of iterations exceeded.")
  return(list(Beta = Beta, Kappa = Kappa, LogL = LogL))
}

#------------------------------------//----------------------------------

formals(LCNRopt) # formals shows the 4 arguments of the function

LeeCarterNR = LCNRopt(dxt   = t(dxt_swiss),
                      ext   = t(ext_swiss), eps = 1e-4, maxiter = 2e3)
Results = list(
  x = ages,
  y = years_swiss,
  beta1  = LeeCarterNR$Beta[, 1],
  beta2  = LeeCarterNR$Beta[, 2],
  kappa2 = LeeCarterNR$Kappa 
)


## We conclude both functions are now working and provide similar results
## We will continue the assignment using the function fit701 from the fitModels script

# Lets now plot the resulting parameter estimates

library(ggplot2)

data_period <- tibble(year = years_swiss, fit = LCfit701$kappa2) ## kappa is time dependent
data_age <- tibble(age = ages_swiss, fit_alpha = LCfit701$beta1, fit_beta = LCfit701$beta2) # both betas are age dependent

## Estimates for beta1 or alpha (age dependent)

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
  geom_line(aes(age, fit_alpha), col = "black") +
  theme_bw() + 
  labs(y = bquote(hat(beta)[x]^"(1)")) 

## Estimates for beta2 or simply beta (age dependent)

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
  geom_line(aes(age, fit_beta), col = "black") +
  theme_bw() +
  labs(y = bquote(hat(beta)[x]^"(2)")) 

## Estimates for kt's (time/year dependent)

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
  geom_line(aes(year, fit), col = "black") +
  theme_bw() + 
  labs(y = bquote(hat(kappa)[t]^"(2)")) 

library(gridExtra)
grid.arrange(g_1, g_2, g_3, ncol = 2, top = "Switzerland - males, 1970 - 2019, Lee Carter, Poisson")


## Plot of Person residuals in a Poisson setting (Residual plot)

grid <- expand.grid(period = years_swiss, age = ages_swiss)
grid$res <- as.vector(LCfit701$epsilon)
names(grid) <- c("Year", "Age", "Residual")

p <- ggplot(grid, aes(x = Year, y = Age)) + geom_tile(aes(fill = Residual)) +
  scale_fill_gradientn(colours =  topo.colors(7)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 15))
p

## It is observable some diagonals, so we see some cohort effects.
# Prof says that diagonal way may indicate that makes sense to include cohort parameters
# in the model, but it usually complicated and lead to unstable results

## -------//------

#Inspecting goodness of fit: observed mx versus fitted mx
# Comparison of observed central death rates and fitted/estimated central death rates
## And we will do it for 4 different ages (25,45,65,75)

age <- 25
rates <- dxt_swiss/ext_swiss ## observed central death rates                                    beta1                          +                 beta2                     *      kt                      
df <- tibble(Year = years_swiss, obs = rates[, age - min(ages_swiss) + 1], fit = exp(LCfit701$beta1[age - min(ages_swiss) + 1] + LCfit701$beta2[age - min(ages_swiss) + 1] * LCfit701$kappa2))

g_25 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + ggtitle("Age 25") + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = bquote(hat(m)[25,][t])) + labs(x = bquote(Year (t)))

age <- 45
rates <- dxt_swiss/ext_swiss
df <- tibble(Year = years_swiss, obs = rates[, age - min(ages_swiss) + 1], fit = exp(LCfit701$beta1[age - min(ages_swiss) + 1] + LCfit701$beta2[age - min(ages_swiss) + 1] * LCfit701$kappa2))

g_45 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + ggtitle("Age45") + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = bquote(hat(m)[45,][t])) + labs(x = bquote(Year (t)))

age <- 65
rates <- dxt_swiss/ext_swiss
df <- tibble(Year = years_swiss, obs = rates[, age - min(ages_swiss) + 1], fit = exp(LCfit701$beta1[age - min(ages_swiss) + 1] + LCfit701$beta2[age - min(ages_swiss) + 1] * LCfit701$kappa2))

g_65 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + ggtitle("Age 65") + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = bquote(hat(m)[65,][t])) + labs(x = bquote(Year (t)))

age <- 75
rates <- dxt_swiss/ext_swiss
df <- tibble(Year = years_swiss, obs = rates[, age - min(ages_swiss) + 1], fit = exp(LCfit701$beta1[age - min(ages_swiss) + 1] + LCfit701$beta2[age - min(ages_swiss) + 1] * LCfit701$kappa2))

g_75 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + ggtitle("Age 75") + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = bquote(hat(m)[75,][t])) + labs(x = bquote(Year (t)))

grid.arrange(g_25, g_45, g_65, g_75, ncol = 2, top = "Switzerland - males, Observed mx vs fitted mx, Lee Carter, Poisson")


## Looking at the graphs we see the model has a very good fit for the ages 45, 65 and 75 
# and an okay fit for the age 25


## Random walk with drift to estimate theta and the variance of the error term
## Use ARIMA (0,1,0)

library(forecast)
time_series = Arima(LCfit701$kappa2, order = c(0,1,0), include.drift = TRUE)

# theta = -2.6522 (negative as it should be)

## Forecast of kt's (more info on my notebook)
## 80,85 and 95 refer to the confidence level
forecast(time_series, level= c(80,85,95))

plot(forecast(time_series, level= c(80,85,95))) # plotting the code line above

# This is how we can quickly calibrate our time series for the Random Walk with drift.

# ALternatively, we can use the function sim2001() of the simModels.R script of the LifeMetrics project
# Now using the simulation function sim2201() (script:simModels.R) of the LifeMetrics

#-----------------------//----------------------  
## Borrow the function sim2001() 

sim2001=function(xx,yy,beta1v,beta2v,kappa2v, mtxLastYear,nsim=100,tmax=20,nyears=0,x0=65, fixStartPoint=FALSE){
  # Model M1
  
  # This simulation does not include parameter uncertainty
  # nsim is the number of sample paths to be generated
  # tmax is the number of years forward that we project
  # nyears is the number of years that we used to fit the random walk
  # x0 is the initial age of a cohort in the final year 
  # that will be projected forward
  
  # inputs xx, yy, kappa1v, kappa2v are key outputs from fit705
  
  # Note the difference between xx and xv defined below
  # xx is the vector of ages used in the original fitting procedure
  #		 and is used as the age range for projecting the 2-dimensional mortality table
  # xv is the vector of ages of the cohort aged x0 in year 1:
  #		 i.e. (x0, x0+1, x0+2, ...)
  
  set.seed(0)
  
  na=length(xx)  # xx=ages used in the fitting procedure
  ny=length(yy)  # yy=calendar years used in the fitting procedure
  
  k2=kappa2v
  if(fixStartPoint) {
    b1 = log(mtxLastYear)
  } else {
    b1=beta1v
  }
  b2=beta2v
  
  d2=diff(k2)	 # differences assumed to be i.i.d.
  m=length(d2)
  # nyears = number of years used to fit the random walk
  # if nyears=0 (default) then all years in data are used
  if(nyears == 0)
  { 
    nyears=m 
  }
  nv=(m-(nyears-1)):m		 # indexes of the last nyears observations
  mu2=mean(d2[nv])		 # mean over the last nyears observations
  v22=mean((d2[nv]-mu2)^2) # MLE for v22
  cc=sqrt(v22)		 	 # standard deviation
  
  n=length(k2)
  
  # tmax is the maximum number of years to follow the cohort
  xv=x0-1+(1:tmax)
  lx=length(xv)
  lxx=length(xx)
  ddv=array(0,c(1,lx))
  # Year 1 correponds to the final year in the sample
  # Hence the year 1 "simulated" value of kappa2 is k2[n]
  dd1=k2[n]
  # ddv is an 1 x lx array with each column being the 
  # values of (k2) over time in a given sample path
  ddv[,1]=dd1
  dda=array(0,c(1,lx,nsim+1))	  # stores the ddv results in a bigger array
  # factor x year x scenario
  tpa=array(0,c(lx,nsim+1))		  # survivor index 
  # year x scenario
  qaa=array(0,c(lxx,tmax,nsim+1)) # q(t,x) mortality rates for all simulations          !!!!!!!!!!!!
  # age x year x scenario
  qa=array(0,c(lxx,tmax))		  # q(t,x) mortality rates for one simulation               !!!!!!!!!!!!
  # age x year
  pa=tpa		 		 		  # one-year survival probabilities, all sims
  # year x scenario
  
  # The Model 1 fitting procedure means that we have no automatic model for
  # mortality outside the age range of xx
  # We use here a crude method for extrapolating beta1 and beta2
  # However, at the end of the simulation loop we trim the 
  # x0 cohort projections (matrix tpa) back to stop at the maximum age xx[na]
  
  # How many extra ages are required to do the full cohort projection?
  extra.ages=max((x0+tmax-1)-xx[na],0)
  # now extrapolate b1 and b2 using a linear projection
  b1ext=c(b1[1:na-1],b1[na]+(b1[na]-b1[1])/(xx[na]-xx[1])*(0:extra.ages) )
  b2ext=c(b2[1:na-1],b2[na]+(b2[na]-b2[1])/(xx[na]-xx[1])*(0:extra.ages) )
  
  # Age indexes for the cohort over the projection
  nv.cohort=(x0-xx[1]+1):(x0-xx[1]+tmax)
  b1cohort=b1ext[nv.cohort]
  b2cohort=b2ext[nv.cohort]
  
  for(e in 1:(nsim+1))
  {   # nsim is the number of sample paths that
    # we want to simulate
    # next loop simulates sample path in years 2 to lx
    for(f in 2:lx)
    {
      if ( e == 1)
        rz = 0
      else {
        rz=rnorm(1)
      }
      
      ddv[,f]=ddv[,f-1]+mu2+cc*rz   ## This line computes the value of kt = kt-1 + theta(mu2) + noise term that we multiply by the std of the normal dist of the error
    }
    dda[,,e]=ddv
    k2v=ddv[1,1:lx]
    if(fixStartPoint) {
      dd1Vector = array(dd1, c(1, length(k2v)))
      mv = exp(b1cohort-b2cohort*(k2v - dd1Vector))
    } else {
      mv=exp(b1cohort+b2cohort*k2v)
    }
    qv=1-exp(-mv)  # vector of mortality rates for this sample path
    # and for the cohort initially aged x0
    pv=(1-qv)	   # convert into survival probabilities
    pa[,e]=pv	   # store in the summary matrix
    tpa[,e]=exp(cumsum(log(pv)))  # vector for the cohort survivor index
    
    # now do the qtx's by age xv[i]
    for(i in 1:lxx)
    {
      if(fixStartPoint) {
        mv0=exp(b1[i] + (b2[i]*(k2v-dd1)))
      } else {
        mv0=exp(b1[i]+b2[i]*k2v)
      }
      qv0=1-exp(-mv0)
      qaa[i,,e]=qv0
    }
  }
  #  Now trim back the tpa matrix
  n1=min(xx[na]-x0+1,tmax)
  pa=pa[1:n1,]
  tpa=tpa[1:n1,]
  
  list(y=yy[n]+(1:tmax)-1,xx=xx,xv=xv[1:n1],mu2=mu2,v22=v22,dda=dda,pa=pa,tpa=tpa,qaa=qaa)
}

#---------------------------------//----------------------------------

# We use 10 000 simulated paths and ask for projections for the next 50 years
# and we use the same number of years used to the Lee Carter calibration for the time series calibration

sim_LC = sim2001(xx = LCfit701$x, yy = LCfit701$y, beta1v = LCfit701$beta1,
                 beta2v = LCfit701$beta2, kappa2v = LCfit701$kappa2, nsim = 10000, tmax =50, 
                 nyears = length(years_swiss))
sim_LC
names(sim_LC)

sim_LC$y # years we are making projections; also includes the last calibrated period 2019
sim_LC$mu2 ## gives the value of theta, again we confirm theta = -2.6522.

sim_LC$dda # gives 10 000 simulated kt's for the next 50 years
dim(sim_LC$dda) # dimensions, we have 10001 because the best estimate is the +1 scenario

sim_LC$dda[ , 1:50, 1] ## the first generated path; Starts with the same kt of 2019; first entry is the best estimate scenario
sim_LC$dda[ , 1:50, 2]## second generated path ; Starts with the same kt of 2019
sim_LC$dda[,,1] ## STORES THE BEST ESTIMATE PATH. why? its basically the 1st path 

LCfit701$kappa2[length(years_swiss)] ## value of the kt in 2019; start of the path

sim_LC$qaa ## gives the one year mortality rates

#--------------------------------//-------------------------------------------------------

# Borrow fan() function

fan=function(yv,mat,pl=0,color){
  # This function plots a fan chart
  # The outer limits are the 5% and 95% quantiles for a given year or age
  # The bands in between are 5% quantiles
  # The darkest band in the middle encompasses the 45% to 55% quantile range
  # Default pl=0 means that the y axis starts at 0 and that the quantity
  # of interest is positive
  # Choose e.g. pl=1 to set the y axis up in a different way
  qv=c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
  ny=length(yv)
  yy=c(yv,yv[ny:1])
  #if(pl == 0){ setaxes(yv[1],yv[ny],0,max(mat)*1.1) }
  #else{ setaxes(yv[1],yv[ny],min(mat),max(mat)) }
  for(i in 1:length(qv)){
    v1=yv*0
    v2=v1
    for(j in 1:ny){
      v1[j]=quantile(mat[j,],qv[i])
      v2[j]=quantile(mat[j,],1-qv[i])
    }
    vv=c(v1,v2[ny:1])
    #plot(c(1961,2054),c(min(k_est),max(k_est)),type="n",main="Kappa_t: estimates + projections (percentiles)",ylab="",xlab="Time")
    # red fan
    if(color=="red"){
      polygon(yy,vv,border=NA,col=rgb(1,(0.5-qv[i])*2,(0.5-qv[i])*2))
    }
    # blue fan
    if(color=="blue"){
      polygon(yy,vv,border=NA,col=rgb((0.5-qv[i])*2,(0.5-qv[i])*2,1))
    }
    # green fan
    if(color=="green"){
      polygon(yy,vv,border=NA,col=rgb((0.5-qv[i])*2,1,(0.5-qv[i])*2))
    }
  }
}

#--------------------------------//-------------------------------------------------------

# time series kappa_t
plot(
  LCfit701$y,
  LCfit701$kappa2,
  pch = 20,
  xlim = c(1970, 2069),
  ylim = range(c(
    range(LCfit701$kappa2), range(sim_LC$dda[, , ])
  )),
  main = bquote(paste("Projection of ", kappa[t])),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(kappa[t])
)

# fan chart
fan(sim_LC$y, sim_LC$dda[, , ], color = "red") ## Include in the report!!!

## In black its the calibrated kt's from 1970 until 2019, and then in red from 2020 until 2068
# These are the projected scenarios for the kt's but of course we can, for each scenario calculate a fan chart for the 
# force of mortality, for death rates qxt,...

# We pick a specific age(here 65) 

age = 65
minage = min(ages_swiss)

plot(
  LCfit701$y,
  exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * LCfit701$kappa2),  ## calculation of the force of mortality for the calibration period(we use LCfit701$kappa2)
  lwd = 2,
  col = "red",
  type = "l",
  ylim = c(0, 0.04),
  main = bquote(paste(
    "Switzerland: Lee-Carter, projection of ", mu[65](t)
  )),
  xlim = c(1970, 2068),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(mu[65](t)))
)
fan(sim_LC$y,
    exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * sim_LC$dda[, , ]),     ## calculation of the force of mortality for the projections(we use sim_LC$dda)
    col = "red")
points(LCfit701$y, rates[, age - minage + 1], col = "black", pch = 20)

## black dots are the force of mortality for the years in the calibration period


## Now we can do the same, but for 1 year death rates, instead of the force of mortality
# We will do for 3 ages:55,65,75

# plot of q_{55,65,75}
ages_sel  = c(55, 65, 75)
nages = length(ages_sel)
color = c("red", "green", "blue")

# y axis is logarithmic
plot(
  c(1970, 2069),
  c(0.002, 0.3),
  type = "n",
  log = "y",
  lwd = 2,
  col = "black",
  main = bquote(paste(
    "Switzerland: projection of ", q[55](t), ", ", q[65](t), ", ", q[75](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[55](t), ", ", q[65](t), ", ", q[75](t)))
)

for (j in 1:nages) {
  lines(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], lwd = 2)
  matlines(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, , 1:1000], col =
             "grey")
}
p = seq(0.05, 0.95, 0.05)

# Construct the corresponding empirical quantiles
for(j in 1:nages) {
  q_int = matrix(data = 0,
                 nrow = 50,
                 ncol = length(p))
  for (i in 1:50)
    q_int[i, ] = quantile(sim_LC$qaa[ages_sel[j] - minage + 1, i, ], p = p)
  # highlight 0.05 and 0.95 quantile
  lines(2019:2068,
        q_int[1:50, 1],
        col = color[j],
        lty = 1,
        lwd = 2)
  lines(2019:2068,
        q_int[1:50, length(p)],
        col = color[j],
        lty = 1,
        lwd = 2)
}

# Alternative: plot of q_{55,65,75} using fan function in simModels.R 

#I prefer this one

plot(
  c(1970, 2058),
  c(0.002, 0.3),
  log = "y",
  main = bquote(paste(
    "Switzerland: projection of ", q[55](t), ", ", q[65](t), ", ", q[75](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[55](t), ", ", q[65](t), ", ", q[75](t))),
  type = "n"
)


for (j in 1:nages) {
  points(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], pch = 20)
  fan(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, ,], color = color[j])
}

##----Kannistö----------------------------------------------------------------------

Switzerland_males_1970_2019_old_age_calib <- Switzerland_males_1970_2019_old_age %>%
  filter(Age >= 80 & Age <= 90)


kan <- lm(data = Switzerland_males_1970_2019_old_age_calib,
          formula = log(mx / (1 - mx)) ~ Age)


phi_1 <- unname(kan$coefficients[1])
phi_2 <- unname(kan$coefficients[2])

mx_old <- (exp(phi_1) * exp(phi_2*c(80:110))) / (1 + exp(phi_1) * exp(phi_2*c(80:110)))

mx_old



og <- Switzerland_males_1970_2019_young %>% 
  filter(Year == 2000) %>% select(c(Age, mx))

kann <- data.frame(c(80:110), mx_old) 

colnames(kann) <- c("Age", "mx")

full <- rbind(og, kann)

ggplot(data = full,
       aes(x = Age, y = log(mx))) +
  geom_line()

##----EPV of 1 unit whole life annuity------------------------------------------------------------








  
