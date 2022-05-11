## Assignment Advanced Life Insurance Mathematics
install.packages("tidyverse")
install.packages("demography")
install.packages("forecast")
library(tidyverse)
library(demography)
library(forecast)



library(dplyr)

## West Germany
#Transforming 110+ into 110

WestGermany_1956_2017 <-
  read_table(file = "./6_2_data/bltper_1x1/DEUTW.bltper_1x1.txt",
             col_names = T,
             skip = 1)

WestGermany_1956_2017["Age"][WestGermany_1956_2017["Age"] =="110+"] <- 110
WestGermany_1956_2017 <- transform(WestGermany_1956_2017,Age = as.numeric(Age))

## Subset from year 1970 until most recent year (2017) (48 years in total)
## Subset from age 0 until 89 because the remaining ages will be closed with Kannisto
WestGermany_1970_2017 = filter(WestGermany_1956_2017, Year >=1970)
WestGermany_1970_2017 = filter(WestGermany_1970_2017, Age <= 89)

years_west <- 1970:max(WestGermany_1956_2017$Year) # May be needed later
ages    <- 0:89 # May be needed later
abc = expand.grid(Year = years_west, Age = ages) # May be needed later
head(abc) # May be needed later

## Plot of the log of the central death rate - we observe an improvement in the central death rates throughout the years

p_west <- ggplot(WestGermany_1970_2017, aes(x = Age, y = log(mx), group = Year)) + 
  geom_line(aes(colour = Year), size = 1, linetype = 1) +
  scale_colour_gradientn(colours = rainbow(10)) +
  scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
  theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") +ggtitle("West Germany - mx evolution")
p_west


abc = expand.grid(Year = years_west, Age = ages)
head(abc)

# Calibration of the Poisson Likelihood (optimization of the Poisson Likelihood 
#with univariate Newton-Raphson steps)

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
    
    # Now apply the constraints
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
#-----------------------//----------------------
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
#-----------------------//----------------------

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
  
#-----------------------//----------------------  

  
ages_west = as.matrix(ages, nrow = 1)
years_west = as.matrix(years_west, nrow = 1)
etx_west = matrix(WestGermany_1970_2017$ex, nrow = 48, byrow= TRUE)
dtx_west = matrix(WestGermany_1970_2017$dx, nrow = 48, byrow= TRUE)

# Estimates for Bx(1), Bx(2), Kt(2)

LCfit701 <-  fit701(ages_west,years_west,etx_west,dtx_west,matrix(1, length(years_west), length(ages)))

names(LCfit701)

LCfit701$beta1 # gives the beta1 for each age (90 ages)
LCfit701$beta2 # gives the beta2 for each age (90 ages)
LCfit701$kappa2 # gives the kt for each year (48 years)

# Lets now plot the resulting parameter estimates

library(ggplot2)

data_period <- tibble(year = years_west, fit = LCfit701$kappa2) ## kappa is time dependent
data_age <- tibble(age = ages_west, fit_alpha = LCfit701$beta1, fit_beta = LCfit701$beta2) # both betas are age dependent

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
  geom_line(aes(age, fit_alpha), col = "black") +
  theme_bw() +
  ggtitle("West Germany - unisex, 1970 - 2017, Lee Carter, Poisson") + 
  labs(y = bquote(hat(beta)[x]^"(1)")) 

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
  geom_line(aes(age, fit_beta), col = "black") +
  theme_bw() + ggtitle("") +
  labs(y = bquote(hat(beta)[x]^"(2)")) 

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
  geom_line(aes(year, fit), col = "black") +
  theme_bw() + ggtitle("") + 
  labs(y = bquote(hat(kappa)[t]^"(2)")) 

library(gridExtra)
grid.arrange(g_1, g_2, g_3, ncol = 2)

## Plot of Person residuals in a Poisson setting (Residual plot)

grid <- expand.grid(period = years_west, age = ages_west)
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

#Inspecting goodness of fit: observed mx versus fitted mx
# Comparison of observed central death rates and fitted/estimated central death rates

age <- 25
rates <- dtx_west/etx_west
df <- tibble(Year = years_west, obs = rates[, age - min(ages_west) + 1], fit = exp(LCfit701$beta1[age - min(ages_west) + 1] + LCfit701$beta2[age - min(ages_west) + 1] * LCfit701$kappa2))

g_25 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.00125, label = "Age 25", size = 10) +
  ggtitle("West Germany - males, 1970 - 2017, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[25,][t])) + labs(x = bquote(Year (t)))

age <- 45
rates <- dtx_west/etx_west
df <- tibble(Year = years_west, obs = rates[, age - min(ages_west) + 1], fit = exp(LCfit701$beta1[age - min(ages_west) + 1] + LCfit701$beta2[age - min(ages_west) + 1] * LCfit701$kappa2))

g_45 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.004, label = "Age 45", size = 10) +
  ggtitle("West Germany - males, 1970 - 2017, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[45,][t])) + labs(x = bquote(Year (t)))

age <- 65
rates <- dtx_west/etx_west
df <- tibble(Year = years_west, obs = rates[, age - min(ages_west) + 1], fit = exp(LCfit701$beta1[age - min(ages_west) + 1] + LCfit701$beta2[age - min(ages_west) + 1] * LCfit701$kappa2))

g_65 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.03, label = "Age 65", size = 10) +
  ggtitle("West Germany - males, 1970 - 2017, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[65,][t])) + labs(x = bquote(Year (t)))

age <- 85
rates <- dtx_west/etx_west
df <- tibble(Year = years_west, obs = rates[, age - min(ages_west) + 1], fit = exp(LCfit701$beta1[age - min(ages_west) + 1] + LCfit701$beta2[age - min(ages_west) + 1] * LCfit701$kappa2))

g_85 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.18, label = "Age 85", size = 10) +
  ggtitle("West Germany - males, 1970 - 2017, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[85,][t])) + labs(x = bquote(Year (t)))

grid.arrange(g_25, g_45, g_65, g_85, ncol = 2)

## Looking at the graphs we see the model has a very good fit for the ages 25,45 and 65, but a very bad 
# fit for the ages of 85.

## Random walk with drift to estimate theta and the variance of the error term
## Use ARIMA (0,1,0)
library(forecast)
time_series = Arima(LCfit701$kappa2, order = c(0,1,0), include.drift = TRUE)

time_series # theta = -2.2004

## Forecast of kt's (more info on my notebook)
## 80,85 and 85 refer to the confidence level
forecast(time_series, level= c(80,85,95))
plot(forecast(time_series, level= c(80,85,95))) # plotting the line above

# Now using the simulation function sim2201() (script:simModels.R) of the LifeMetrics
## Here we will also see longevity charts

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
  qaa=array(0,c(lxx,tmax,nsim+1)) # q(t,x) mortality rates for all simulations
  # age x year x scenario
  qa=array(0,c(lxx,tmax))		  # q(t,x) mortality rates for one simulation
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
      
      ddv[,f]=ddv[,f-1]+mu2+cc*rz
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
# Now using the simulation function sim2201() (script:simModels.R) of the LifeMetrics
## Here we will also see longevity charts
# We use 10 000 simulated paths and ask for projections for the next 50 years
# and we use the same number of years used to the Lee Carter calibration for the time series calibration

sim_LC = sim2001(xx = LCfit701$x, yy = LCfit701$y, beta1v = LCfit701$beta1,
      beta2v = LCfit701$beta2, kappa2v = LCfit701$kappa2, nsim = 10000, tmax =50, 
      nyears = length(years_west))
sim_LC
names(sim_LC)

sim_LC$y # years we are making projections; also includes the last calibrated period 2017

dim(sim_LC$dda) # dimensions, we have 10001 because the best estimate is the +1 scenario

sim_LC$dda[ , 1:50, 1] ## the first generated path; Starts with the same kt of 2017
sim_LC$dda[ , 1:50, 2]## second generated path ; Starts with the same kt of 2017
sim_LC$dda[,,1] ## STORES THE BEST ESTIMATE PATH

LCfit701$kappa2[length(years_west)] ## value of the kt in 2017; start of the path

## fan chart 


##see slides 41 onwards

#-----------------------//---------------------------------------------





