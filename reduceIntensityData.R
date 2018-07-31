library(sp)
library(ggmap)
library(raster)

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

# set wd
setwd("/home/mhaas/PhD/Routines/IPE_comparison")

#Read in data set
data <- read.csv("DESERVE_Intensity.csv")

#set all data without depth to default_depth
data$depth[which(is.na(data$depth))] <- default_depth

#create spatial inputs
sites <- SpatialPoints(coords=data.frame(data$lon, data$lat),proj4str=crsys)
epicentres <- SpatialPoints(coords=data.frame(data$epi_lon, data$epi_lat),proj4str=crsys)

#calculate epicentral distance for comparison
R_epi <- calcRepi(sites,epicentres)
data$R_epi <- R_epi


# ##################################################
# # Create fake sites for given R_hyp from (lon,lat)=(0,0) (for comparison)
# ##################################################
# createSites <- function(R_hyp){
#   #just go from (0,0) to east
#   lat <- rep(0,length(R_hyp))
#   lon <- R_hyp/(2*pi*6371/360)
#   sites <- SpatialPoints(coords=data.frame(lon,lat),proj4str=crsys)
#   return(sites)
# }
# 
# ##########################
# # Check probability with bindi ipe
# ##########################
# 
# fsites<-createSites(c(150))
# epi <- SpatialPoints(coords=data.frame(rep(0,1),rep(0,1)) ,proj4str=crsys)
# est <- BPOAMZ2011(7fsites,epi,type='epi')

########################
# Reduction of dataset
#########################

#1. remove rows observations > 200km
data <- data[which(R_epi < 200),]

#2. remove C.E. events (23 same as catalogue)
data <- data[which(data$year>=23),]

#3.remove observations > 7.5 in 150km
idxs1 <- which(data$R_epi>150)
idxs2 <- which(data$MeanInt >= 7.5)
idxs <- intersect(idxs1,idxs2)

data<-data[-idxs,]

#4. restrict to intensity V+
data <- data[which(data$MeanInt >= 5),]
# 
# data$check <- rep(0,length(R_epi))
# data$check[idxs] <- 1

write.csv(data,'reduced_intensity_data.csv')
