library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)

source("kaklamanos_ngaR_SourceCode/Misc.R")
source("kaklamanos_ngaR_SourceCode/resorce.R")
source("kaklamanos_ngaR_SourceCode/NGA.R")
source("kaklamanos_ngaR_SourceCode/BA08.R")
source("kaklamanos_ngaR_SourceCode/AS08.R")
source("kaklamanos_ngaR_SourceCode/CB08.R")
source("kaklamanos_ngaR_SourceCode/CY08.R")

Vs30<-800  
epsilon<-0
rake <- NA
U<-0
SS<-0 
NS<-1 
RS<-0
AB11<-1
Z2.5 <- 3500
Ztor <- 0 
dip <- NA
rake <- 90
Rrup <- NA
Rx<- NA
azimuth<- NA
W<- NA
Fhw<-0

mag<-c(4,5,6,7)
nmag<-length(mag)
metri<-c(1,10,20,50,100,200)
nmetri<-length(metri)
ncl<- nmetri*nmag

modelNGAW1 <- array(0, c(1,ncl), list(NULL,c(1:ncl)))
modelNGAW1 <- list()

# create NGA08 models
ii<-0
for (ik in 1:nmetri)
{
  for (ij in 1:nmag)
  {
    print(mag[[ij]])
    ii<-ii+1
    modelNGAW1[[ii]] <- Sa.nga(M = mag[[ij]], dip = dip, rake = rake, Fhw = Fhw, Ztor = Ztor, W = W, 
                             Rjb = metri[[ik]], Rrup = Rrup, Rx = Rx, azimuth = azimuth, 
                             Vs30 = Vs30, epsilon = 0, T = c(4,1,0.1,0.02))
  }
}
#create resorce BINDI model
modelRESO <- array(0, c(1,ncl), list(NULL,c(1:ncl)))
modelRESO <- list()
T<-c(4,1,0.1,0.02)
ii<-0
for (ik in 1:nmetri)
{
  for (ij in 1:nmag)
  {
    print(mag[[ij]])
    ii<-ii+1
    modelRESO[[ii]] <- Sa.bi(M = mag[[ij]],Rjb = metri[[ik]],Vs30 = Vs30, epsilon = 0, T = c(4,1,0.1,0.02),SS,NS,RS)
  }
}

# GO FOR SAMMON'S MAPS
out.sam <- list()
p1 <- list()
for (jj in 1:4)  #periods
{
     pro<-matrix(0,5,ncl)
     for (ik in 1:ncl){
          pro[1,ik]<-modelNGAW1[[ik]]$Y50.as[jj]
          pro[2,ik]<-modelNGAW1[[ik]]$Y50.ba[jj]
          pro[3,ik]<-modelNGAW1[[ik]]$Y50.cb[jj]
          pro[4,ik]<-modelNGAW1[[ik]]$Y50.cy[jj]
          pro[5,ik]<-modelRESO[[ik]][jj]
     }
     proT <- cbind(pro,c(1:5))
     
     dist.data <- dist(proT[,1:ncl])
     randstart <- matrix(runif(nrow(proT[,1:ncl])*2),nrow(proT[,1:ncl]),2)
     out.sam[[jj]] <- sammon(dist.data,y=randstart,tol=1e-8,niter=500)
     dum <- array(NA, c(dim(out.sam[[jj]]$points)[1], 3), list(NULL,c("x","y","GMPE")))
     modsom <-data.frame(dum)
     modsom$x <- out.sam[[jj]]$points[,1]
     modsom$y <- out.sam[[jj]]$points[,2]
     modsom$class <- c("AS08","BA08","CB08","CY08","BI14")
     
     p1[[jj]]<-ggplot()+
       geom_point(data= modsom,aes(x,y,color=factor(class)),pch=19,cex=6)+
       scale_x_continuous(name="X",limits=c(-0.5,0.5))+scale_y_continuous(name="Y",limits=c(-0.5,0.5))+
       ggtitle(paste("Sammon's plot for T=",T[jj],"s"))
}#periods

#PLOT
p1[[1]]
p1[[2]]
p1[[3]]
p1[[4]]
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]])# ,ncol = 2, nrow=2, main ="")
