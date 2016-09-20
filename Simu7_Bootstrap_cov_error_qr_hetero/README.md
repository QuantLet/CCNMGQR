
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **Bootstrap_cov_error_qr_hetero** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : Bootstrap_cov_error_qr_hetero

Published in : Confidence corridors for nonparametric multivariate generalized quantile regression

Description : 'The main file which performs the Monte Carlo simulation for the coverage ratio of
the bootstrap multivariate confidence band for nonparametric kernel quantile regression. The data
generating model is a heterogeneous model.'

Keywords : kernel, nonparametric, multivariate, quantile, regression, confidence-bands, bootstrap

See also : kernel, lcrq

Author : Shih-Kang Chao, Katharina Proksch, Holger Dette and Wolfgang Haerdle

Submitted : 19/09/2016 by Lining Yu

```


### R Code:
```r
# clear all variables
rm(list = ls(all = TRUE))
graphics.off()
# set the working directory
#setwd("C:/...")
libraries = c("np", "quantreg", "VGAM", "rgl", "misc3d", "matrixStats", "MASS")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
source("lcrq.r")
source("kernel.r")

#########################     General setting       ##############################

f<-function(x1,x2){sin(2*pi*x1)+x2}
S <- matrix(c(1, -0.3, -0.3, 1), nrow = 2) 
bb=seq(0.1,0.9,length=20)
xx<-as.matrix(expand.grid(bb,bb))
Rep <- 500
tau50<-0.5
tau20<-0.2
tau80<-0.8
q50<-qnorm(tau50)
q20<-qnorm(tau20)
q80<-qnorm(tau80)
sig_0 <- 0.2 # Please change here for the other model variance: 0.5 0.7
nn <- 100 # Please change here for the other sample sizes: 300, 500

################     Hetergeneity   sigma = f(x)     #################################

sig <- function(x1,x2){sig_0+0.8*x1*(1-x1)*x2*(1-x2)}
f0_hom_tau50<-f(xx[,1],xx[,2])+q50*sig(xx[,1],xx[,2])
f0_hom_tau20<-f(xx[,1],xx[,2])+q20*sig(xx[,1],xx[,2])
f0_hom_tau80<-f(xx[,1],xx[,2])+q80*sig(xx[,1],xx[,2])

error_hom<-matrix(0,nrow=3,ncol=3)
area.cc <- numeric(0)
temp.area<-c(0,0,0)
########         n = 50              ################################################ 

for(k in 1:Rep){

X <- mvrnorm(nn, mu = c(0,0), Sigma = S)
x.biunif <- pnorm(X) 
y.biunif<-f(x.biunif[,1],x.biunif[,2])+rnorm(nn,mean=0,sd=sig(xx[,1],xx[,2]))

bdwh<-npcdensbw(xdat=x.biunif,ydat=y.biunif,bwmethod="normal-reference",ckertype="epanechnikov",ckerorder=2)
bdwh

try(qr.ufit50<-lcrq.boot(x.biunif,y.biunif,bwth=bdwh,d=2,tau=tau50,xx=xx))
try(qr.ufit20<-lcrq.boot(x.biunif,y.biunif,bwth=bdwh,d=2,tau=tau20,xx=xx))
try(qr.ufit80<-lcrq.boot(x.biunif,y.biunif,bwth=bdwh,d=2,tau=tau80,xx=xx))

try(if(length(which(qr.ufit50$lband - f0_hom_tau50>0 | qr.ufit50$hband - f0_hom_tau50<0))!=0){error_hom[1,1]<-error_hom[1,1]+1})
try(if(length(which(qr.ufit20$lband - f0_hom_tau20>0 | qr.ufit20$hband - f0_hom_tau20<0))!=0){error_hom[1,2]<-error_hom[1,2]+1})
try(if(length(which(qr.ufit80$lband - f0_hom_tau80>0 | qr.ufit80$hband - f0_hom_tau80<0))!=0){error_hom[1,3]<-error_hom[1,3]+1})

temp.area<-c(0,0,0)
try(temp.area[1]<-sum(qr.ufit50$hband-qr.ufit50$lband))
try(temp.area[2]<-sum(qr.ufit20$hband-qr.ufit20$lband))
try(temp.area[3]<-sum(qr.ufit80$hband-qr.ufit80$lband))
try(area.cc<-rbind(area.cc,temp.area))

print(k)
}

error_hom/Rep
```
