#Sammon Map for IPEs
rm(list=ls(all.names=TRUE))
library(MASS)
library(ggplot2)
library(sp)
library(raster)
library(kohonen)

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

##################################################
# Create fake sites for given R_epi from (lon,lat)=(0,0) (for comparison)
##################################################
 createSites <- function(R_epi){
   #just go from (0,0) to east
   lat <- rep(0,length(R_epi))
   lon <- R_epi/(2*pi*6371/360)
   sites <- SpatialPoints(coords=data.frame(lon,lat),proj4str=crsys)
   return(sites)
 }

#######################
# IPEs
#######################
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

########################
# Sammons map
########################

#1) discrete parameter space for model comparison
mags <- seq(5,7.5,.1)
repi <- seq(5,200,5)
sites <- createSites(repi)
epis <- SpatialPoints(coords=data.frame(rep(0,length(sites)),rep(0,length(sites))) ,proj4str=crsys)

#2) get estimates for all models
IPEs <- data.frame(SSG2009,AWW2012,SSG2010,CL2002,LBB2014,BPOAMZ2011,B2006,DR2005,PAGDL2008,GVTB2010)
magnitudes <- matrix(0,length(IPEs),length(mags)*length(repi))
repis <- matrix(0,length(IPEs),length(mags)*length(repi))
means <- matrix(0,length(IPEs),length(mags)*length(repi))
sigmas <- matrix(0,length(IPEs),length(mags)*length(repi))
for (i in seq(length(IPEs))){
  started <- 0
  #go through all magnitudes
  for (mag in mags){
    estimate <- IPEs[[i]](rep(mag,length(sites)),sites,epis)
    if (started==0){
      tmp_magnitudes <- rep(mag,length(sites))
      tmp_repi <- repi
      tmp_means <- estimate$value.int
      tmp_sigmas <- estimate$value.sigmas
      started<-1
    }else{
      tmp_magnitudes <- c(tmp_magnitudes,rep(mag,length(sites)))
      tmp_repi <- c(repi,tmp_repi)
      tmp_means <- c(tmp_means,estimate$value.int)
      tmp_sigmas <- c(tmp_sigmas,estimate$value.sigmas)
    }
  }
  magnitudes[i,]<-tmp_magnitudes
  repis[i,]<-tmp_repi
  means[i,]<- tmp_means
  sigmas[i,]<- tmp_sigmas
}

sam<-sammon(dist(means),tol=1e-8,niter=500)
sam$class <- c('SSG2009','AWW2012','SSG2010','CL2002','LBB2014','BPOAMZ2011','B2006','DR2005','PAGDL2008','GVTB2010')
x <- sam$points[,1]
y <- sam$points[,2]
df <- data.frame(x,y,ipe=sam$class)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#CCCCCC")
h <- ggplot()+geom_point(data=df,aes(x,y,color=ipe))+scale_colour_manual(values=cbPalette)
h

#clustering
cluster <- kmeans(cbind(means,sigmas),5)
ss <- cluster$withinss
cl <- cluster$cluster
ss_out<-rep(99,length(cl))
for (i in seq(length(cl))){
  ss_out[i]<-ss[cl[i]]
}
df <- data.frame(x,y,ipe=sam$class,cluster=cl,ss=ss_out)

write.csv(df,file='sammon.map')


#hc <- hclust(dist(cbind(x,y)),method="ward.D")


# ###########################
# # Kohonen Map
# ###########################
# #reorder data as Matrix with columns Mw,R,mean,sigma
# data <- cbind(as.vector(magnitudes),as.vector(repis),as.vector(means))#,as.vector(sigmas))
# 
# # Create a training data set (rows are samples, columns are variables
# # Here I am selecting a subset of my variables available in "data"
# set.seed(42)
# data_train <- data[sample(seq(length(data[,1])),1000),]
# 
# # Change the data frame with training data to a matrix
# # Also center and scale all variables to give them equal importance during
# # the SOM training process. 
# data_train_matrix <- as.matrix(scale(data_train))
# 
# # Create the SOM Grid - you generally have to specify the size of the 
# # training grid prior to training the SOM. Hexagonal and Circular 
# # topologies are possible
# som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")
# 
# # Finally, train the SOM, options for the number of iterations,
# # the learning rates, and the neighbourhood are available
# som_model <- som(data_train_matrix,grid=som_grid,rlen=100,alpha=c(0.05,0.01),keep.data = TRUE,n.hood='circular' )
# 
# ################
# # Plot models on SOM
# ###############
# model1 <- cbind(magnitudes[1,],repis[1,],means[1,])
# #model1 <- cbind(means[1,],sigmas[1,])
# #model1 <- cbind(means[1,])
# model1.sc <- scale(model1)
# model1.som <- predict(som_model,trainX=model1[,1:2],trainY=model1[,3])
# plot(som_model,type='mapping',classif=model1.som)
# 
# data(wines)
# set.seed(7)
# 
# training <- sample(nrow(wines), 120)
# Xtraining <- scale(wines[training, ])
# Xtest <- scale(wines[-training, ],center = attr(Xtraining, "scaled:center"),scale = attr(Xtraining, "scaled:scale"))
# 
# som.wines <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))
# 
# som.prediction <- predict(som.wines, newdata = Xtest,trainX = Xtraining,trainY = factor(wine.classes[training]))
# table(wine.classes[-training], som.prediction$prediction)
# 



# #########################
# # Check dist calculation
# #########################
# first <- c(1,2,3,4,5,6,7,8,9)
# second <- c(2,2,2,4,5,6,7,8,8)
# third <- c(9,8,7,6,5,4,3,2,1)
# m <- rbind(first,second,third)
# dist(m,diag=T,upper=T)