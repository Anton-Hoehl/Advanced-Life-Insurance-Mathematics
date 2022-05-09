## Finland males

install.packages("tidyverse")
install.packages("demography")
install.packages("forecast")
library(tidyverse)
library(demography)
library(forecast)

#Transforming 110+ into 110
Finland_males_1878_2020["Age"][Finland_males_1878_2020["Age"] =="110+"] <- 110
Finland_males_1878_2020 <- transform(Finland_males_1878_2020,Age = as.numeric(Age))

## Subset from year 1970 until most recent year (2019) (50 years in total)
## Subset from age 0 until 89 because the remaining ages will be closed with Kannisto
Finland_males_1970_2019 = filter(Finland_males_1878_2020, Year >=1970 & Year <=2019)
Finland_males_1970_2019 = filter(Finland_males_1970_2019, Age <= 89)

years_fin <- 1970:max(Finland_males_1970_2019$Year) # May be needed later
ages    <- 0:89 # May be needed later
abc = expand.grid(Year = years_west, Age = ages) # May be needed later
head(abc) # May be needed later

## Plot of the log of the central death rate - we observe an improvement in the central death rates throughout the years

p_fin <- ggplot(Finland_males_1970_2019, aes(x = Age, y = log(mx), group = Year)) + 
  geom_line(aes(colour = Year), size = 1, linetype = 1) +
  scale_colour_gradientn(colours = rainbow(10)) +
  scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
  theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") +ggtitle("Finland Male population - mx evolution")
p_fin

# Calibration of the Poisson Likelihood (optimization of the Poisson Likelihood 
#with univariate Newton-Raphson steps)

source('fitModels.R')

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


ages_fin = as.matrix(ages, nrow = 1)
years_fin = as.matrix(years_fin, nrow = 1)
etx_fin = matrix(Finland_males_1970_2019$ex, nrow = 50, byrow= TRUE)
dtx_fin = matrix(Finland_males_1970_2019$dx, nrow = 50, byrow= TRUE)

# Estimates for Bx(1), Bx(2), Kt(2)

LCfit701 <-  fit701(ages_fin,years_fin,etx_fin,dtx_fin,matrix(1, length(years_fin), length(ages_fin)))

names(LCfit701)

LCfit701$beta1 # gives the beta1 for each age (90 ages)
LCfit701$beta2 # gives the beta2 for each age (90 ages)
LCfit701$kappa2 # gives the kt for each year (48 years)






