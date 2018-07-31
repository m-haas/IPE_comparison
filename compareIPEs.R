rm(list=ls(all.names=TRUE))
library(ggplot2)
library(ggmap)
library(reshape)
library(raster)
library(arm)

###########################################################################
# Compare different IPE Performance to historic data
# 1. Implemented IPEs
# 2. Read in data
#    a) pre 1900 event data (Sbeinati et al 2005)
#    b) 1927 event data (Zohar and Marco 2012)
#    c) 1995 event data (Al-Tarazi 2000) 
# 3. Bayes Inference: Goodness of fit IPE -> Data a,b,c
# 4. Regression from Data
###########################################################################
#global constants
#crs
crsys <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#depth
default_depth <- 10

###########################################################################
# Functions
###########################################################################
calcRhyp <- function(sites,epicentres,depth){
  # Calculates hypocentral distances for given sites and hypocentres
  # input:  1) sites (set of spatial points)
  #         2) epicentres (set of spatial points)
  #         3) depths [km]
  # output: 1) hypocentral distances vector [km]
  R_epi <- c()
  for (i in seq_along(sites)){
    R_epi <- append(R_epi,spDists(sites[i],epicentres[i],longlat=TRUE))
  }
  return(sqrt(R_epi^2+depth^2))
}
calcRepi <- function(sites,epicentres){
  # Calculates epicentral distances for given sites and epicentre
  # input:  1) sites (set of spatial points)
  #         2) epicentres (set of spatial point)
  # output: 1) epicentral distances vector [km]
  R_epi <- c()
  for (i in seq_along(sites)){
    R_epi <- append(R_epi,spDists(sites[i],epicentres[i],longlat=TRUE))
  }
  return(R_epi)
}

###########################################################################
# Available IPEs
###########################################################################
SSG2009 <- function(mags,sites,epicentres,depths=default_depth,type='epi',R_JB=c(0)){
  # Intensity Prediction equation according to
  # Bibtex: Sorensen2009
  # 2 different types of distance measure epicentral/joyner-boore
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  #           4) depth = hypocentral depth [km]
  #           5) type = distance type ['epi','JB']
  #           6) Joyner&Boore distance for each event (if type == 'JB')
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity[0=no,1=yes]
  
  if (type=='JB'){
    a <- 0.376
    b <- 5.913
    c <- -2.656
    d <- -0.0020
    sigma <- 0.672
    #distance measure is Joyner and Boore distance
    R <- R_JB
  }else{
    a <- 0.793
    b <- 3.417
    c <- -2.157
    d <- -0.0065
    sigma <- 0.742
    #distance measure is epicentral
    R <- calcRepi(sites,epicentres)
  }
  #valid ranges
  mrange <- c(5.9,7.4)
  rrange <- c(0,350)
  irange <- c(5,10)
  
  #calculate intensities
  int <- a*mags + b + c*log10(sqrt((R^2+depths^2)/depths^2)) + d*(sqrt(R^2+depths^2)-depths)
  sigmas <- rep(sigma,length(int))
  
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <- 0
  valid[which(int < irange[1]| int > irange[2])] <- 0
  
  #return the intensity estimates, the sigmas and a flag (1 valid, 0 invalid) if the value is extrapolated, i.e. invalid
  return(data.frame(int,sigmas,valid))
}

AWW2012 <- function(mags,sites,epicentres,depths=default_depth){
  # Calculates the intensities resulting from an event at given site locations
  # implemented from Allen, Wald & Worden 2012
  # only the hypocentral distance model was implemented (not the rupture distance model)
  # input: 1) mag = Moment magnitude (Scalar)
  #        2) sites = locations to calculate intensity at (Set of spatial points)
  #        3) epicentre = epicentre (Spatial point)
  #        4) depth = hypocentral depth [km]
  # output:1) dataframe containing magnitude,intensity[MMI],sigma[MMI],epicentral distance [km]
  
  #Parameters for IPE
  c0 <-  2.085
  c1 <-  1.428
  c2 <-  -1.402
  c4 <-  0.078
  m1 <-  -0.209
  m2 <-  2.042
  #Parameters for distance dependent uncertainty
  s1 <- 0.82
  s2 <- 0.37
  s3 <- 22.9
  
  mrange <- c(5.0,7.9)
  rrange <- c(0,300)
  irange <- c(1,12)
  #Calculate hypocentral distances
  R <- calcRhyp(sites,epicentres,depths)
  
  #Calculate attenuation
  #Rm
  Rm <- m1 + m2 * exp(mags-5)
  #close site effect
  indSwitchON <- which(R>50)
  switchON <- matrix(0,length(R),1)
  switchON[indSwitchON] <- 1
  #return intensities,sigma and validity (ln(dist))
  int <- c0+c1*mags+c2*log(sqrt(R^2+Rm^2))+switchON*c4*log(R/50)
  sigmas <- s1 + s2/(1+(R/s3)^2)
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

SSG2010 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Sorensen2010
  # only the one using epicentral distance and Standard regression was implemented
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag
  
  c <-  1.556
  e <- -0.428
  a <-  5.518
  b <- -0.0020
  h <-  15.550
  sigma <- 0.972
  
  mrange <- c(6.3,7.0)
  rrange <- c(0,300)
  irange <- c(3,11)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  #estimate intensities  
  int <- c*mags + e - a*log10(sqrt((R^2+h^2)/h^2)) - b*(sqrt(R^2+h^2)-h)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

CL2002 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Chandler2001
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <- -0.8919
  b <-  1.4798
  c <-  0.1311
  d <- -0.0364
  e <-  0.0193
  f <-  0.0085
  sigma <- 0.5 #not precisely specified (p783)
  #valid ranges
  mrange <- c(3.3,8.0)
  rrange <- c(0,300)
  irange <- c(4,10)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  R0 <- 0.5*10^(0.74*mags-3.55)
  #distance dependent effects
  switch1 <- rep(0,length(R))
  switch2 <- switch1
  switch1[which(R > 45 && R <= 75)] <- 1
  switch2[which(R > 75)] <- 1
  #estimate intensities  (ln!)
  int <- a + b*mags - c*log(sqrt((R+R0)/R0)) + d*R + switch1*e*(R-45) + switch2*f*(R-75)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

LBB2014 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: LeGoff2014
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <- -1.9438
  b <-  4.1
  c <- -9.5763
  sigma <- 0.63
  
  #valid ranges
  mrange <- c(4.4,6.2)
  rrange <- c(0,1000)
  irange <- c(1,11)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  
  #estimate intensities  (ln!)
  int <- a*log(R) + b*mags + c
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

BPOAMZ2011 <- function(mags,sites,epicentres,depths=default_depth,type='hyp'){
  # Intensity Prediction equation according to
  # Bibtex: Bindi2011
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  #           4) type = ['hyp','epi','fix'] hypocentral,epicentral,fixed depth (15km)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  #convert Mw to Ms (since IPE derived with MLH) Surface wave magnitude conversion according to Scordilis2006
  mags[which(mags<6.2)]<-(mags[which(mags<6.2)]-2.07)/0.67
  mags[which(mags>=6.2)]<-(mags[which(mags>=6.2)]-0.08)/0.99
  
  if(type=='hyp'){
    a1 <-  1.071
    a2 <-  1.003
    a3 <-  2.621
    a4 <-  5.567*10^(-4)
    sigma <- 0.710
    R <- calcRhyp(sites,epicentres,depths)
    int = a1*mags + a2 - a3*log10(R/10) - a4*(R-10) 
  }else if(type=='epi'){
    a1 <-  0.898
    a2 <-  1.215
    a3 <-  1.809
    a4 <-  3.447*10^(-3)
    sigma <- 0.737
    R <- calcRepi(sites,epicentres)
    int = a1*mags + a2 - a3*log10(sqrt((R^2+depths^2)/depths^2)) - a4*(sqrt(R^2+depths^2)-depths)
  }else{
    a1 <-  1.049
    a2 <-  0.686
    a3 <-  2.706
    a4 <-  1.811*10^(-4)
    sigma <- 0.689
    R <- calcRepi(sites,epicentres)
    int = a1*mags + a2 - a3*log10(sqrt((R^2+depths^2)/depths^2)) - a4*(sqrt(R^2+depths^2)-depths)
  }
  
  #valid ranges
  mrange <- c(4.6,8.3)
  rrange <- c(0,600)
  irange <- c(4,10)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

B2006 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Bakun2006
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <-  0.44
  b <-  1.70
  c <- -0.0048
  d <- -2.73
  h <- 10
  sigma <- 0.58
  
  #valid ranges
  mrange <- c(4.6,7.3)
  rrange <- c(0,500)
  irange <- c(3,8)
  
  #distance measure is hypocentral distance with fixed depth h
  R <- calcRhyp(sites,epicentres,h)
  
  #estimate intensities  
  int <- a + b*mags + c*R +d*log10(R)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

DR2005 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Dowrick2005
  # only the main seismic region one (exists also for different focal mechanisms e.g. strike-slip,reverse etc)
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a1 <-  4.40
  a2 <-  1.26
  a3 <- -3.67
  a4 <-  0.012
  a5 <-  0.409
  d  <- 11.78
  sigma <- 0.43
  
  #valid ranges
  mrange <- c(4.6,8.2)
  rrange <- c(0,500)
  irange <- c(3,11)
  
  #distance measure is rupture distance approximated by hypocentral distance
  R <- calcRhyp(sites,epicentres,depths)
  
  #switch for crustal event(we consider only crustal)
  dc <- 1
  #estimate intensities
  int <- a1 + a2*mags + a3*log10((R^3+d^3)^(1/3.)) + a4*depths + a5*dc
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

PAGDL2008 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Pasolini2008
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <-  0.0086
  b <-  1.037
  #epicentral int from GOR and using only measured Mw not derived Mw
  c <- -4.446 #typo in publication according to Gasperini2010
  d <-  2.210 #typo in publication according to Gasperini2010
  
  h <- 3.91
  sigma <- 0.69
  
  #valid ranges
  mrange <- c(4.4,7.4)
  rrange <- c(0,200)
  irange <- c(4,11)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  D <- sqrt(R^2+h^2)
  
  #epicentral intensity
  IE <- c + d*mags 
  
  #estimate intensities (ln!!)
  int <- IE - a*(D-h) - b*(log(D)-log(h))
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

GVTB2010 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Gasperini2010
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <-  0.0009
  b <-  1.172
  #epicentral int
  c <- -5.368 #typo in publication according to Gasperini2010
  d <-  2.364
  
  h <- 4.49
  sigma <- 0.74
  
  #valid ranges
  mrange <- c(4.4,7.4)
  rrange <- c(0,1000)
  irange <- c(4,11)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  D <- sqrt(R^2+h^2)
  
  #epicentral intensity
  IE <- c + d*mags 
  
  #estimate intensities
  int <- IE - a*(D-h) - b*(log(D)-log(h))
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  return(data.frame(int,sigmas,valid))
}
##################################################
#
##################################################
#own IPE based on regression
ownIPE <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: 
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <- 0.126027
  b <- 1.239118
  c <- 0.001066
  d <- -1.902045
  eta <- 0.8063
  sigma <- 0.6956
# # no-random effect but using median of distance binned 
#   a <- 0.92531
#   b <- 0.98655
#   c <- -0.01409
#   d <- 0.72880
#   sigma <- 1.101861
  
  R <- calcRhyp(sites,epicentres,depths)
  R2 <- R-10
  R <- R/10
  int <- a + b*mags + c*R2 + d*log10(R)
  
  #valid ranges
  mrange <- c(4.0,7.5)
  rrange <- c(0,500)
  irange <- c(3,8)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

IPEs <- data.frame(SSG2009,AWW2012,SSG2010,CL2002,LBB2014,BPOAMZ2011,B2006,DR2005,PAGDL2008,GVTB2010)


##################################################
# Observation statistics
##################################################
# Plots histogram of provided observed intensity for different events
# destinguished by EventID and providing year of occurrence
StatPlot <- function(filename,EventID,year,obs){
  nr_obs <- length(obs)
  nr_events <- length(unique(EventID))
  period <- c(min(year),max(year))
  text <- paste(nr_obs,'observations for',nr_events,'events between',period[1],'and',period[2],sep=' ')
  df <- data.frame(obs)
  h <- ggplot(df,aes(obs))
  h + geom_histogram(binwidth = 1) + 
      scale_x_discrete(name = 'Intensity', breaks=seq(round(min(obs)),round(max(obs)))) + 
      ylab('Frequency') + 
      ggtitle(text)
  ggsave(filename)
}

# # Plots Intensity and Distance pairs
# IntDistPlot <- function(filename,dist_epi,int){
#   df <- data.frame(dist_epi,int)
#   h <- ggplot(df,aes(x=dist_epi,y=int))+
#     geom_point() + 
#     scale_x_log10(name = 'Epicentral distance') + 
#     scale_y_continuous(name='Intensity',limits=c(3,11),breaks=c(3,4,5,6,7,8,9,10,11),labels=c('III','IV','V','VI','VII','VIII','IX','X','XI'))
#     ggtitle(text)
#   ggsave(filename)
# }

# # Plot epicentres and Intensity Observations
# EpiObsPlot <- function(filename,epicentres,sites){
#   bb <- bbox(sites)
#   epicentres <- data.frame(epicentres)
#   sites <- data.frame(sites)
#   colnames(epicentres)<-c('lon','lat')
#   epicentres <- unique(epicentres)
#   colnames(sites)<-c('lon','lat')
#   map <- ggmap(get_map(location = bb)) + geom_point(data=sites,aes(lon,lat),shape=2,show_guides=TRUE) + 
#         geom_point(data=epicentres,colour='red',shape=8,show_guides=TRUE)
#   ggsave(filename)
# }
##################################################
# Goodness of fit
##################################################
#1) residuals
MRP <- function(IPE_name,observed,estimated,distance){
  #Mean residual plot
  #Generates a plot which shows the residual (obs-predicted)
  #the values were ordered by epicentral distance
  #permutation<-order(distance)
  #observed <- observed[permutation]
  #estimated <- estimated[permutation]
  residuals <- observed-estimated
  Repi <- distance
  #Repi <- distance[permutation]
  df <- data.frame(Repi,residuals)
  p <- ggplot(df, aes(Repi,residuals))
  p + geom_point(aes(Repi,residuals)) +
      ggtitle(IPE_name) +
      ylab("Residual") +
      xlab("Epicentral distance")
  ggsave(paste(IPE_name,'residuals.eps',sep='_'))
  
  #return sum of square residuals
  return(sqrt(mean(residuals^2)))
}

#2) EDR Kale and Akkar 2013
EDR <- function(observed,estimated,sigmas,dd=0.01,x=3){
  # Adopted algortihm from original MATLAB script
  # moving from lognormal to normal distibution for macroseismic intensity
  # Euclidean-Distance Based Ranking (EDR) Method
  # This script is used for performing EDR Method.
  # EDR Method is developed by Ozkan Kale and Sinan Akkar.
  # ozkankale@gmail.com; sakkar@metu.edu.tr
  # 05.04.2012
  # Earthquake Engineering Research Center
  # Department of Civil Engineering
  # Middle East Technical University
  # 06800, Ankara, Turkey
  
  # Format of Input File ("InputData.txt"):
  # 
  # Intensity Measure File -> 1st column: ln(Observed)
  #                           2nd column: ln(Estimated)
  #                           3rd column: Corresponding sigma values
  # 
  # 
  # Input values
  # dd = 0.01;          % infinitesmall bandwidth
  # x = 3;              % the multiplier of sigma
  # 
  # Intensity Measure File
  # SpecQuantFile = fopen('InputData.txt','r');
  # data = fscanf(SpecQuantFile,'%g\n',[3 inf]);
  # data = data';
  a <- observed         #ln of observed ground motion intensity
  dmin <- dd/2          #minimum discrete distance (euclidean distance model-observed)
  
  # Calculation of the Kappa value
  #mean of observations
  mua <- mean(a)
  
  #estimates of the attenuation model and their mean
  Y <- estimated
  muY <- mean(Y)
  
  # Assuming a log-log linear relation
  # Calculate the coefficients of Y=B1*A+B0 (Fitted line)
  x1 <- (a - mua*(Y - muY))
  y1 <- (a - mua^2)
  B1 <- sum(x1)/sum(y1)
  B0 <- muY - B1*mua
  
  Yfit <- B0 + B1*a
  Yc <- Y - (Yfit - a)      # corrected estimations (Equation 10)
  
  DEorg <- sqrt(sum((a-Y)^2))          # Equation 9.b
  DEcor <- sqrt(sum((a-Yc)^2))         # Equation 9.c
  
  Kappa <- DEorg/DEcor        # Equation 9.a
  
    # Calculation of the EDR value
  mu_Y <- estimated
  s_Y <- sigmas #not log see definition of residuals, i.e. (log(x_obs)-log(x_est))/sigma_est
  
  mu_D <- a - mu_Y            # Equation 3.a
  s_D <- s_Y                  # Equation 3.b
  
  # Selection of an appropriate dmax value (Equation 8)
  d1c <- abs(a - (mu_Y - x*s_Y))
  d1cmax <- max(d1c)
  d2c <- abs(a - (mu_Y + x*s_Y))
  d2cmax <- max(d2c)
  dcmax <- ceiling(max(d1cmax, d2cmax))
  
  dmax <- dcmax - dd/2         # selected maximum discrete distance
  nd <- length(seq(dmin,dmax,dd))   # number of discrete increment
  
  MDE <- 0
  for (j in seq(1,nd)){
    di <- dmin+(j-1)*dd
    d <- di
    d1i <- di-dd/2
    d1 <- d1i
    d2i <- di+dd/2
    d2 <- d2i
    
    # Calculations given in Equation 5
    P1 <- pnorm((d1-mu_D),s_D) - pnorm((-d1-mu_D),s_D)
    P2 <- pnorm((d2-mu_D),s_D) - pnorm((-d2-mu_D),s_D)
    
    # Calculations given in Equation 6
    P <- P2-P1
    MDEi <- P*d
    MDE <- MDE + MDEi 
  }
  
  # Modified Euclidean Distance normalized by N (The first component of EDR index in terms of MDE)
  MDE_norm <- sqrt(1/length(observed)*sum(MDE^2))
  # Square root of Kappa (The second component of EDR index in terms of Kappa)
  Kappa_sq <- sqrt(Kappa)
  # Resultant EDR value (Equation 11)
  EDR <- sqrt(Kappa*1/length(observed)*sum(MDE^2))
  
  # Return the results
  return(data.frame(MDE_norm,Kappa_sq,EDR))
}

#3) Log-likelihood - Scherbaum,Delavaud and Riggelsen 2009
LLH <- function(observed,estimated,sigmas){
  # average log-likelihood described by Scherbaum et al. 2009
  # assumption that returned value of the gmpe is mean of a 
  # normally distributed groundmotion pdf at this distance and 
  # sigmas is the standard deviation of the pdf 
  
  # normally distributed pdf at each location 
  g <- dnorm(observed,mean=estimated,sd=sigmas)
  # return average sample log-likelihood in units of bits (log2)
  #return (1/length(observed)*sum(log2(g)))
  return (1/length(observed)*sum(log2(g)))
}

KLD <- function(LLH1,LLH2){
  # Relative Kullback-Leibler (KL) distance between two models
  # approximated by the difference between their average log-likelihood
  
  # return relative KL-distance
  # i.e. relative information loss for model 2 compared to model 1
  return (LLH2-LLH1)
}

KLW <- function(LLHs){
  #Returns Kullback-Leibler weights for the given log-likelihoods
  # in units of bits (log2)
  avg <- sum(2^(LLHs))
  weights <- 2^(LLHs)/avg
  
  return(weights)
}

##################################################
# Create fake sites for given R_hyp from (lon,lat)=(0,0) (for comparison)
##################################################
createSites <- function(R_hyp){
  #just go from (0,0) to east
  lat <- rep(0,length(R_hyp))
  lon <- R_hyp/(2*pi*6371/360)
  sites <- SpatialPoints(coords=data.frame(lon,lat),proj4str=crsys)
  return(sites)
}

#################
# Comparison
#################
# set wd
setwd("/home/mhaas/PhD/Routines/IPE_comparison")

#Read in data set
#data <- read.csv("all_intensity_data.csv")
data <- read.csv("reduced_intensity_data.csv")

#only 1927,1995
#data <- data[which(data$year > 1900),]

#set all data without depth to default_depth
data$depth[which(is.na(data$depth))] <- default_depth

#create spatial inputs
sites <- SpatialPoints(coords=data.frame(data$lon, data$lat),proj4str=crsys)
epicentres <- SpatialPoints(coords=data.frame(data$epi_lon, data$epi_lat),proj4str=crsys)

#calculate epicentral distance for comparison
R_epi <- calcRepi(sites,epicentres)

#average observed intensity (renew for each)
observed <- data$MeanInt

#create histogram of observed intensity
StatPlot('intensity_histogram.eps',data$EventID,data$year,observed)
#create Int-Dist Plot
#IntDistPlot('int_dist.eps')

#create model estimates
rankStats <- data.frame("IPE" = character(1),
                        RMS = numeric(1),
                        LLH = numeric(1),
                        MDE_norm = numeric(1),
                        Kappa_sq = numeric(1),
                        EDR = numeric(1),
                        stringsAsFactors=FALSE)
for (i in seq_along(IPEs)){
  estimate <- IPEs[[i]](data$Mw,sites,epicentres,depths=data$depth)
  #use only these which are valid (according to author)
  idx <- which(estimate$value.valid==1)
  #execute if at least one is valid
  if (length(idx)>=1){
    int <- estimate$value.int[idx]
    obs <- observed[idx]
    Repi <- R_epi[idx]
    sig <- estimate$value.sigmas[idx]
  
    #generate mean residuals plot for IPE and return RMS
    RMS <- MRP(colnames(IPEs[i]),obs,int,Repi)
    #calculate MDE
    EDR_res <- EDR(obs,int,sig,dd=0.01,x=3)
    #calculate log-likelihood
    LLH_res <- LLH(obs,int,sig)
    #append to data.frame
    rankStats <- rbind(rankStats,c(colnames(IPEs)[i],RMS,LLH_res,EDR_res$MDE_norm,EDR_res$Kappa_sq,EDR_res$EDR))
  }
}
#remove first fake row (added to don't raise warnings)
rankStats <- rankStats[-c(1),]
#get a ranking of the IPEs RMS,abs(LLH) and EDR are equally weighted
points <- rank(as.numeric(rankStats['RMS'][[1]])) +
          rank(abs(as.numeric(rankStats['LLH'][[1]]))) +
          rank(as.numeric(rankStats['EDR'][[1]]))
#rank the IPEs according to their points (fewest=1st etc.)
ranking <- rankStats['IPE'][order(points),]

# #choose first n IPEs and calculate weights
# n <- 3
# llhs <- as.numeric(rankStats['LLH'][[1]][order(points)])[1:n]
# w <- KLW(llhs)

#create latex ready table with stats and ranking
out <- data.frame(seq(1:length(IPEs)),rankStats['IPE'][order(points),],as.numeric(rankStats['RMS'][order(points),]),as.numeric(rankStats['LLH'][order(points),]),as.numeric(rankStats['EDR'][order(points),]),as.numeric(c(0.5,0.25,0.25,0,0,0,0,0,0,0)))
out[,3:6]<-round(out[-c(1,2)],2)
colnames(out) <- c('Rank','IPE','RMSR','LLH','EDR','Weigth')
write.table(out,file='ranks.latex',round(digits=2),sep='&',quote=FALSE,row.names=FALSE,col.names=TRUE)
re <- scan('ranks.latex',what='character')
out <- c()
for (row in re){
  out<-rbind(out,paste(row,'\\',sep=''))
}
write.table(out,file='ranks.latex',quote=FALSE,row.names=FALSE,col.names=FALSE)


################################
# Regression model
################################

library(raster)
library(reshape)

data$R_hyp <- calcRhyp(sites,epicentres,data$depth)
# bin the available data by magnitude
mag_bin <- 0.5
mmin <- floor(min(data$Mw)*1/mag_bin)*mag_bin
mmax <- ceiling(max(data$Mw)*1/mag_bin)*mag_bin
bins <- seq(mmin,mmax-mag_bin,mag_bin)
#get indices of events for each bin and store bin in the data.frame
data$bin <- seq_along(data$Mw)
for (bin in bins){
  idxs <- which(data$Mw >= bin & data$Mw < (bin+mag_bin))
  data$bin[idxs] <- bin
}

#function to plot int,dist for specific mag bin compare to some estimate function
plotBinPairs <- function(bin,binwidth,data,fcts,fctnames,fname){
  idxs <- which(data$bin >= bin)
  #create fake sites for comparison
  x <- seq(1,1000,1)
  site <- createSites(x)  
  epi <- SpatialPoints(coords=data.frame(rep(0,1000),rep(0,1000)) ,proj4str=crsys)
  #z <- data$depth[idxs]
  # fix depth for estimate
  z <- 10
  x <- calcRhyp(site,epi,z)
  #remove distance duplicates
  i2keep <- !duplicated(x)
  x <- x[i2keep]
  site <- site[i2keep,]
  epi <- epi[i2keep,]
  #z <- z[i2keep] 
  z <- rep(z,length(epi))
  #order by distance
  iOrder <- order(x)
  x <- x[iOrder]
  site <- site[iOrder,]
  epi <- epi[iOrder,]
  z <- z[iOrder]
  #get estimates and plot as lines
  mags <- rep(bin+binwidth/2,length(epi))
  #create df
  df <- data.frame()
  for (i in seq_along(fcts)){
    fct <- fcts[[i]]
    fctname <- fctnames[i]
    y <- fct(mags,site,epi,z)
    mean <- y$int
    sd <- y$sigmas
    name <- rep(fctname,length(mean))
    tdf <- data.frame(x,mean,sd,name)
    df <- rbind(df,tdf) 
  }
  #write IPE estimates to gmt csv
  write.table(df,paste(fname,'_gmt_est.csv',sep=''),sep=',')
  #write observation to gmt2.csv
  write.table(data[idxs,],paste(fname,'_gmt_obs.csv',sep=''),sep=',')
  # center <- bin+binwidth/2
  # text <- paste('Mw_sim:',as.character(center),'\n Mw_obs:',as.character(center),'+/-',as.character(binwidth/2))
  p <- ggplot()+ 
       #standard deviation as ribbon
       geom_ribbon(data=df,aes(x=x,ymax=mean+sd,ymin=mean-sd,fill=name,colour=name),linetype=2,alpha=0.5)+
       #mean estimates
       geom_line(data=df,aes(x=x,y=mean,colour=name),size=1)+
       #annotate("text", x = 250, y = 11, label = text,size=3)+
       #mean observed
       geom_point(data=data[idxs,],aes(R_hyp,MeanInt),show_guide=TRUE)+
       scale_x_continuous(limits=c(10,300))+
       scale_y_continuous(limits=c(3,11),breaks=c(3,4,5,6,7,8,9,10,11),labels=c('III','IV','V','VI','VII','VIII','IX','X','XI'))+
       ylab('Intensity')+
       xlab('Hypocentral Distance [km]')
       
  #save plot
  ggsave(paste(fname,'_Mw7_25z10km.eps',sep=''))
}

#Hough&Avni
HA2009 <- function(mags,sites,epicentres,depths=default_depth){
  # Intensity Prediction equation according to
  # Bibtex: Hough2009
  # input:    1) mags = Moment magnitudes
  #           2) sites = locations to calculate intensity at (Set of spatial points)
  #           3) epicentres = epicentres (Set of spatial points)
  # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
  
  a <-  -0.64
  b <-  1.70
  c <- -0.00448
  d <- -1.67
  sigma <- 0.0
  
  #valid ranges
  mrange <- c(4.6,7.3)
  rrange <- c(0,500)
  irange <- c(3,8)
  
  #distance measure is epicentral distance
  R <- calcRepi(sites,epicentres)
  
  #estimate intensities  
  int <- a + b*mags + c*R +d*log10(R)
  
  # standard deviation
  sigmas <- rep(sigma,length(int))
  # check validity: mrange,rrange,irange
  valid <- rep(1,length(mags))
  valid[which(mags < mrange[1]| mags > mrange[2])] <- 0
  valid[which(R < rrange[1]| R > rrange[2])] <-0
  valid[which(int < irange[1]| int > irange[2])] <-0
  
  return(data.frame(int,sigmas,valid))
}

#plot IPEs against the observations for each cluster
plotBinPairs(7.,mag_bin,data,c(SSG2009,B2006,DR2005,SSG2010),c('SSG2009','B2006','DR2005','SSG2010'),'cluster1')
plotBinPairs(7.,mag_bin,data,c(CL2002),c('CL2002'),'cluster3')
plotBinPairs(7.,mag_bin,data,c(PAGDL2008,GVTB2010),c('PAGDL2008','GVTB2010'),'cluster2')
plotBinPairs(7.,mag_bin,data,c(AWW2012,BPOAMZ2011),c('AWW2012','BPOAMZ2011'),'cluster4')
plotBinPairs(7.,mag_bin,data,c(LBB2014),c('LBB2014'),'cluster5')

###########################
# actual regression model
###########################
# restrict to data R_hyp <= 300 km and sample 2/3 for regression and 1/3 for testing
#data2 <- data[which(data$R_hyp <= 300),]

# #bin data by distance
# bin <- 10
# dists <- seq(0,300,bin)
# dists2 <- c()
# medInt <- c()
# mags <- seq(4,8,.1)
# mags2 <- c()
# len <- c()
# depths <- c()
# #get median Int for each distance bin
# for (m in mags){
#   idxs1 <- which(data2$magnitude == m)
#   for (d in dists){
#     idxs2 <- which(data2$R_hyp >= d & data2$R_hyp < d+bin)
#     idxs <- intersect(idxs1,idxs2)
#     if(length(idxs)!=0){
#       len <- c(len,length(idxs))
#       mags2 <- c(mags2,m)
#       dists2 <- c(dists2,d)
#       medInt <- c(medInt,median(data2$MeanInt[idxs]))
#       depths <- c(depths,data2$depth[idxs])
#     }
#   }
# }
# #regression
# #Bindi style
# Mw <- mags2
# R <- dists2/10
# R2 <- dists2-10
# I <- medInt
# df <- data.frame(Mw,R,R2,I)
# reg <- glm(I ~ Mw + R2 + log10(R),data=df)

# #Dowrick style
# Mw <- mags2
# R <- dists2
# I <- medInt
# df <- data.frame(I,Mw,R)
# 
# reg <- glm(I ~ Mw + R+ log10(R),data=df)


#sampling deactivated

#nr <- round(length(data[,1])/3)
#samples <- sample(seq(1,length(data[,1])),nr)
#test_data <- data[samples,]
#data2 <- data[-samples,]
#differentiate by Rhyp
# data1$dc <- rep(1,length(data2[,1]))
# data2$dc[which(data2$R_hyp > 50)] <- 2
#  
# R2 <- data2$R_hyp-10
# reg1 <- glm(MeanInt ~ magnitude + R2 + log10(R_hyp/10),family='gaussian',data2)
# reg2 <- lmer(MeanInt ~ magnitude + R2 + log10(R_hyp/10) + (1|EventID),data2)
# reg3 <- lmer(MeanInt ~ magnitude + R2 + log10(R_hyp/10) + (1|EventID) + (1|dc),data2)
# reg4 <- lmer(MeanInt ~ magnitude + R2 + log10(R_hyp/10) + magnitude*R_hyp + (1|EventID) + (1|dc),data2)
#  
# #AWW style
# reg5 <- lmer(MeanInt ~ magnitude + log(sqrt(R_hyp^2+(1+exp(magnitude-5))^2))+(1|EventID),data2)
# reg6 <- lmer(MeanInt ~ magnitude + log(sqrt(R_hyp^2+(1+exp(magnitude-5))^2))+(1|EventID),data2)
# # 
# #get performance values
# IPE <- function(mags,sites,epicentres,depths){
# # #   a1 <- 3.33243
# # #   a2 <- 0.77013
# # #   a3 <- -0.00402
# # #   a4 <- -1.32231
# #   a1 <- -0.011595
# #   a2 <- 1.292655
# #   a3 <- -0.004334
# #   a4 <- -1.3509
# #   sigma <- 1.0729
#   a1 <- 0.779517
#   a2 <- 1.185166
#   a3 <- -0.004016
#   a4 <- -1.585783
#     
#   R <- calcRhyp(sites,epicentres,depths)
#   int = a1 + a2*mags + a3*log10(R/10) + a4*(R-10)
#   # standard deviation
#   sigmas <- rep(sigma,length(int))
#   # check validity
#   valid <- rep(1,length(int))
#   return(data.frame(int,sigmas,valid))
# }
# plotBinPairs(7.,.5,data,c(IPE),c('IPE'))

# # get starting values
# data2$Rmin10 <- data2$R_hyp-10
# data2$Rdiv10 <- data2$R_hyp/10
# reg <- glm(MeanInt ~ magnitude + Rmin10 + log10(Rdiv10),data=data2)
# theta <- unname(coef(reg))
# start <- list(theta)
# reg <- lmer(MeanInt ~ magnitude + Rmin10 + log10(Rdiv10)+(1|EventID),data=data2,start=start)
# 
# #confidence intervals by bootstrapping
# reg.boot = confint.merMod(reg,c(1,2,3,4,5,6),level=0.95,method="boot",nsim=1000,FUN=NULL,boot.type="norm")
# #F-tests

# #catterpillar plot random effect
# ## re = object of class ranef.mer
# ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
#   f <- function(x) {
#     pv   <- attr(x, "postVar")
#     cols <- 1:(dim(pv)[1])
#     se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
#     ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
#     pDf  <- data.frame(y=unlist(x)[ord],
#                        ci=1.96*se[ord],
#                        nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
#                        ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
#                        ind=gl(ncol(x), nrow(x), labels=names(x)))
#     
#     if(QQ) {  ## normal QQ-plot
#       p <- ggplot(pDf, aes(nQQ, y))
#       p <- p + facet_wrap(~ ind, scales="free")
#       p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
#     } else {  ## caterpillar dotplot
#       p <- ggplot(pDf, aes(ID, y)) + coord_flip()
#       if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
#         p <- p + facet_wrap(~ ind)
#       } else {           ## different scales for random effects
#         p <- p + facet_grid(ind ~ ., scales="free_y")
#       }
#       p <- p + xlab("Levels") + ylab("Random effects")
#     }
#     
#     p <- p + theme(legend.position="none")
#     p <- p + geom_hline(yintercept=0)
#     p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
#     p <- p + geom_point(aes(size=1.2), colour="blue") 
#     return(p)
#   }
#   
#   lapply(re, f)
# }
# ggCaterpillar(ranef(reg, condVar=TRUE))

# #bootstrapping analysis
# FUN <- function(fit) {
#   return(fixef(fit))
# }
# reg.boot = bootMer(reg,FUN, nsim = 10)

# # manifold regression using 80% data 20% testing
# set.seed = 42
# sampleNr = round(length(data2$magnitude)*0.8)
# #vectors to store
# as <- c()
# bs <- c()
# cs <- c()
# ds <- c()
# sigs <- c()
# RMSs <- c()
# LLHs <- c()
# EDRs <- c()
# for (i in seq(1,10)){
#   set.seed = 42+i
#   #split observations (idxs) in derivation and test observations
#   idxs <- seq_along(data2$magnitude)
#   idxs_der <- sample(idxs,sampleNr)
#   idxs_tst <- setdiff(idxs,idxs_der)
#   der <- data2[idxs_der,]
#   tst <- data2[idxs_tst,]
#   #run regression on sample using start values
#   reg <- lmer(MeanInt ~ magnitude + Rmin10 + log10(Rdiv10)+(1|EventID),data=der,start=start)
#   #define IPE
#   IPE <- function(mags,sites,epicentres,depths=default_depth){
#       # Intensity Prediction equation according to
#       # Bibtex: 
#       # input:    1) mags = Moment magnitudes
#       #           2) sites = locations to calculate intensity at (Set of spatial points)
#       #           3) epicentres = epicentres (Set of spatial points)
#       # output:1) dataframe containing intensity[MMI],sigma[MMI],validity flag  
#       
#       a <- fixef(reg)[1]
#       b <- fixef(reg)[2]
#       c <- fixef(reg)[3]
#       d <- fixef(reg)[4]
#       sigma <- sigma(reg)
#       # # no-random effect but using median of distance binned 
#       #   a <- 0.92531
#       #   b <- 0.98655
#       #   c <- -0.01409
#       #   d <- 0.72880
#       #   sigma <- 1.101861
#       
#       R <- calcRhyp(sites,epicentres,depths)
#       R2 <- R-10
#       R <- R/10
#       int <- a + b*mags + c*R2 + d*log10(R)
#       sigmas <- rep(sigma,length(int))
# 
#       return(data.frame(int,sigmas))
#     }
#   #run performance test using test data
#   sites <- SpatialPoints(coords=data.frame(tst$lon, tst$lat),proj4str=crsys)
#   epicentres <- SpatialPoints(coords=data.frame(tst$epi_lon, tst$epi_lat),proj4str=crsys)
#   estimate <- IPE(tst$magnitude,sites,epicentres,depths=tst$depth)
#   est <- estimate$int
#   sig <- estimate$sigmas
# 
#   #generate mean residuals plot for IPE and return RMS
#   RMS <- MRP("IPE",tst$MeanInt,est,tst$R_epi)
#   #calculate MDE
#   EDR_res <- EDR(tst$MeanInt,est,sig,dd=0.01,x=3)
#   #calculate log-likelihood
#   LLH_res <- LLH(tst$MeanInt,est,sig)
#   #store estimated parameters and the performance
#   as <- c(as,fixef(reg)[1])
#   bs <- c(bs,fixef(reg)[2])
#   cs <- c(cs,fixef(reg)[3])
#   ds <- c(ds,fixef(reg)[4])
#   sigs <- c(sigs,sigma(reg))
#   RMSs <- c(RMSs,RMS)
#   LLHs <- c(LLHs,LLH_res)
#   EDRs <- c(EDRs,EDR_res$EDR)
#   #set start values a new
#   theta <- unname(fixef(reg))
#   start <- list(theta)
# }
# #combine results
# df <- data.frame(as,bs,cs,ds,sigs,RMSs,LLHs,EDRs)