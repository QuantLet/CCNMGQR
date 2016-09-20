# clear all variables
rm(list = ls(all = TRUE))
graphics.off()
# set the working directory
#setwd("C:/...")
libraries = c("CBPS", "boot", "cubature", "np", "quantreg", "VGAM", "rgl", "matrixStats")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
source("kernel.r")
source("lcrq.r")

############################ Retrieve data ########################################################### 
data(LaLonde)
lalonde<-LaLonde[2491:3212,]  # The males

la.treat<-lalonde[which(lalonde$treat==1),] # treatment group
la.contl<-lalonde[which(lalonde$treat==0),] # control group

x.tr<-cbind(la.treat$age,la.treat$educ)
x.co<-cbind(la.contl$age,la.contl$educ)
######################################################################################################

########################### Data arrangement #########################################################
y.tr.r75<-la.treat$re75/1000   # in thousand dollar
y.co.r75<-la.contl$re75/1000
y.tr.r78<-la.treat$re78/1000
y.co.r78<-la.contl$re78/1000
y.tr<- y.tr.r78-y.tr.r75 # Dependent variable of treatment group: income growth in three years 78-75
y.co<- y.co.r78-y.co.r75 # Dependent variable of control group: income growth in three years 78-75

xx1<-seq(from=19,to=31,by=1)
xx2<-seq(from=7,to=13,by=1)
xx<-as.matrix(expand.grid(xx1,xx2))
########################## Choose bandwidths ##########################################################
bdwh.tr<-npcdensbw(xdat=x.tr,ydat=y.tr,ckertype="epanechnikov",ckerorder=2) # Cross-validated bandwidth obtained from function npcdensbw in package "np"
bdwh.tr
bdwh.tr.80<-bdwh.tr
bdwh.tr.90<-bdwh.tr
bdwh.tr.80$xbw<-bdwh.tr$xbw*1.3
bdwh.tr.90$xbw<-bdwh.tr$xbw*1.7
bdwh.co<-npcdensbw(xdat=x.co,ydat=y.co,ckertype="epanechnikov",ckerorder=2) # ckertype="epanechnikov"
bdwh.co
########################## Quantile regression for treatment and control group for different quantile levels ####################################################
system.time(
boot50.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr,d=2,tau=0.5,xx=xx,B=10000)
)
system.time(
boot50.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.5,xx=xx,B=10000)
)
system.time(
boot30.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.3,xx=xx,B=10000)
)
system.time(
boot70.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.7,xx=xx,B=10000)
)
system.time(
boot20.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.2,xx=xx,B=10000)
)
system.time(
boot80.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.8,xx=xx,B=10000)
)
system.time(
boot90.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.9,xx=xx,B=10000)
)
system.time(
boot10.co<-lcrq.boot(x=x.co,y=y.co,bwth=bdwh.co,d=2,tau=0.1,xx=xx,B=10000)
)
#############
system.time(
boot30.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.3,xx=xx,B=10000)
)
system.time(
boot70.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.7,xx=xx,B=10000)
)
system.time(
boot20.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.2,xx=xx,B=10000)
)
system.time(
boot80.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.8,xx=xx,B=10000)
)
system.time(
boot90.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.90,d=2,tau=0.9,xx=xx,B=10000)
)
system.time(
boot10.tr<-lcrq.boot(x=x.tr,y=y.tr,bwth=bdwh.tr.90,d=2,tau=0.1,xx=xx,B=10000)
)
################# Present results by pictures #################################################
#dev.off()
persp3d(xx1,xx2,boot10.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-50,30))  #
persp3d(xx1,xx2,boot10.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot10.tr$fv,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-50,30))
persp3d(xx1,xx2,boot10.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot10.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot10.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot10.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot20.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-45,30))  #
persp3d(xx1,xx2,boot20.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot20.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot20.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot20.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot20.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot30.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-12,12))  #
persp3d(xx1,xx2,boot30.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot30.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot30.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot30.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot30.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot50.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-25,45))  #
persp3d(xx1,xx2,boot50.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot50.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot50.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot50.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot50.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot70.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-25,45))  #
persp3d(xx1,xx2,boot70.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot70.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot70.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot70.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot70.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot80.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-25,45))  #
persp3d(xx1,xx2,boot80.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot80.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot80.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot80.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot80.pdf",fmt="pdf")

#dev.off()
persp3d(xx1,xx2,boot90.tr$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-75,100))  #
persp3d(xx1,xx2,boot90.tr$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray46", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot90.co$hband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot90.co$lband,alpha = 1,lwd=0.6,color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_boot90.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_boot90.pdf",fmt="pdf")
####################################################################################################################################################################################################################################################
# The sample size of treatment group is smaller, so we use local linear kernel estimator to estimate the quantile functions of treatment group again, and compare them with the confidence bands for control group
llfit10.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.90,d=2,tau=0.1,xx=xx)
llfit20.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.2,xx=xx)
llfit30.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.3,xx=xx)
llfit50.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr,d=2,tau=0.5,xx=xx)
llfit70.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.7,xx=xx)
llfit80.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.80,d=2,tau=0.8,xx=xx)
llfit90.tr<-llrq(x=x.tr,y=y.tr,bwth=bdwh.tr.90,d=2,tau=0.9,xx=xx)


persp3d(xx1,xx2,boot10.tr$fv,alpha=1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-50,30))
persp3d(xx1,xx2,boot10.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot10.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr10.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr10.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit20.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-12,12))
persp3d(xx1,xx2,boot20.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot20.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr20.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr20.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit30.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-12,12))
persp3d(xx1,xx2,boot30.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot30.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr30.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr30.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit50.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-5,20))
persp3d(xx1,xx2,boot50.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot50.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr50.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr50.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit70.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-5,20))
persp3d(xx1,xx2,boot70.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot70.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr70.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr70.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit80.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-5,20))
persp3d(xx1,xx2,boot80.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot80.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr80.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr80.pdf",fmt="pdf")

persp3d(xx1,xx2,llfit90.tr$fv,alpha = 1,lwd=0.6,color="gray46",xlab = "Age", ylab = "Schooling in years",zlab = "Earnings in 78-75", box=FALSE,axes=FALSE,xlim = c(19,29),ylim =c(7,13),zlim=c(-25,45))
persp3d(xx1,xx2,boot90.co$hband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
persp3d(xx1,xx2,boot90.co$lband,alpha = 1,lwd=0.6,back="line",front="line",color="gray95", xlab = "", ylab = "",zlab = "",add=TRUE)
axes3d(edges=c('x--','y--','z--'))

   #rgl.snapshot("lalonde_fvtr90.png", fmt="png", top=TRUE )
   rgl.postscript("lalonde_fvtr90.pdf",fmt="pdf")

