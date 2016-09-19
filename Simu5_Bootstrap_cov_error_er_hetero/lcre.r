# Depend: np, matrixStats, VGAM, ks


lcre<-function(x, y, H, d, tau = 0.5, xx)  # H: sequence of bandwidth, xx: a matrix m*d, x: matrix n*d, y: n*1 vector
{
    fv <- numeric(0)
    # dv <- xx
    n <- length(y)
    for (i in 1:length(xx[,1])) {
        ll<-matrix(1,d,n)
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- epanech.prod(z[k,]/H,d)
        }
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    list(xx = xx, fv = fv)
}
##################################################################################
loss.phi<-function(x,tau){
phii<-numeric(0)
phii<-(2*abs(x)*(as.numeric(x<0)-tau))^2  # "squared" loss
return(phii)
}
###############################################################################
phi <- function(x,tau){
phii<-numeric(0)
phii<-2*abs(x)*(as.numeric(x<0)-tau)
return(phii)
}
###############################################################################
lcre.conf<-function(x, y, bwth, d, tau = 0.5, kern="quartic", alpha=0.05,xx){
if(kern=="epanech"){
K_2 <- (3/5)^d # Product kernel Epanechnikov
SIG <- diag(d)*(3/2)*(3/5)^(d-1)  # Product kernel Epanechnikov
ker <- epanech.prod
}
if(kern=="quartic"){
K_2 <- (5/7)^d # Product kernel quartic
SIG <- diag(d)*(15/7)^(5/7)^(d-1)  # Product kernel quartic
ker <- quartic.prod
}

fv <- numeric(0)
fv.i <- numeric(0)  # for estimating the conditional varepsilon^2(x)
sigma.x <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.01)   # Undersmoothing
    ll<-matrix(1,d,n)
    for (i in 1:length(xx[,1])) {
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    i <- 1
    k <- 1
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv.i[i] <- r@coefficients[1]
    }
   eps.hat <- y-fv.i
   phi.i<-loss.phi(eps.hat,tau=tau)
    i <- 1
    k <- 1
     for (i in 1:length(xx[,1])){
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- lm(phi.i ~ 1, weights = wx)
        sigma.x[i] <- r$coefficients[1]
    }
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=d)     
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   Fe_x.objt<-npcdist(txdat=x,tydat=y,bwmethod="normal-reference",exdat=xx,eydat=fv)
   Fe_x<-Fe_x.objt$condist
   xx.mat<-as.matrix(xx)
   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
   kap <- abs(log(prod(H)^(1/d), base = n))
   eta <- log(vol.D)+ kap*d*log(n)
   c_alpha <- log(2)-log(abs(log(1-alpha))) 
   H_2 <- vol.D*prod(H/H[1])*(2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))
   band <- 1/sqrt(n*prod(H))*sqrt(sigma.x*sqrt(K_2)/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*(d_n+c_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, sigma.x=sigma.x, Fe_x=Fe_x,eps.hat=eps.hat)
}
###################################################################################
lcre.boot<-function(x, y, bwth, d, tau = 0.5, kern="quartic", B=500, alpha=0.05,xx){
if(kern=="epanech"){
K_2 <- (3/5)^d # Product kernel Epanechnikov
SIG <- diag(d)*(3/2)*(3/5)^(d-1)  # Product kernel Epanechnikov
ker <- epanech.prod
}
if(kern=="quartic"){
K_2 <- (5/7)^d # Product kernel quartic
SIG <- diag(d)*(15/7)^(5/7)^(d-1)  # Product kernel quartic
ker <- quartic.prod
}

fv <- numeric(0)
fv.i <- numeric(0)  # for estimating the conditional varepsilon^2(x)
sigma.x <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.01)   # Undersmoothing
    ll<-matrix(1,d,n)
    W.x <- matrix(0,nrow=length(xx[,1]),ncol=n)
    for (i in 1:length(xx[,1])) {
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        W.x[i,] <- wx
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    i <- 1
    k <- 1
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv.i[i] <- r@coefficients[1]
    }
   eps.hat <- y-fv.i
   phi.i<-loss.phi(eps.hat,tau=tau)
    i <- 1
    k <- 1
     for (i in 1:length(xx[,1])){
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- lm(phi.i ~ 1, weights = wx)
        sigma.x[i] <- r$coefficients[1]
    }
   
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=d)     
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   Fy_x.objt<-npcdist(txdat=x,tydat=y,bwmethod="normal-reference",exdat=xx,eydat=fv)   
   Fe_x.objt<-npcdist(txdat=x,tydat=eps.hat,bwmethod="normal-reference",exdat=xx,eydat=seq(0,length=length(xx[,1]),by=0))   
   Fe_x<-Fe_x.objt$condist
   Fy_x<-Fy_x.objt$condist   
#   xx.mat<-as.matrix(xx)
#   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
#   kap <- abs(log(prod(H)^(1/d), base = n))
#   eta <- log(vol.D)+ kap*d*log(n)
#   H_2 <- vol.D*prod(H/H[1])*(2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
#   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))
#####################################################################################################################
   RX<-diag(1/sqrt(sigma.x*fx),nrow=length(fx),ncol=length(fx))  #*((Fe_x*(2*tau-1)-tau)/(Fy_x*(2*tau-1)-tau))       Diagonal matrix with 1/sqrt(f(x))being the diagonal line
   BOOT <- cbind(x,eps.hat)  # epsilon is at the last column of this matrix
   sample.boot <- numeric(0) # bootstrap sample
   HH<-H
   HH[d+1]<-Fe_x.objt$ybw
   stoch.dev <- matrix(0,nrow=length(xx[,1]),ncol=B)  # storage for bootstrap deviation
   max.dev.boot <- as.vector(seq(0,length=B,by=0))
   W.x.b <- matrix(0,nrow=length(xx[,1]),ncol=n)
   i <- 1
   k <- 1
   l <- 1
   # system.time(
   for(b in 1:B){
   bb<-sample(c(1:n),size=n,replace=TRUE)   
   er<-rnorm(n*(d+1))
   E<-matrix(er,nrow=n,ncol=d+1)
   boot.data<-as.matrix(BOOT[bb,])+ E%*%diag(HH,nrow=d+1,ncol=d+1)
  # xx.hat<-cbind(xx,seq(0,0,length=length(xx[,1])))

   x.boot <- boot.data[,1:d]
   sample.boot <- boot.data[,d+1]
   psi.sample.boot <- as.matrix(phi(sample.boot,tau))
   W.x<-matrix(0,nrow=length(xx[,1]),ncol=n)
   X.boot<-matrix(t(x.boot),nrow=n*length(xx[,1]),ncol=d,byrow=TRUE)
   XX<-xx[rep(seq_len(nrow(xx)),each=n),]   
   Z<-X.boot-XX
   Eff<- which(abs(t(t(Z)/H))[,1]<1 & abs(t(t(Z)/H))[,2]<1)    ## Only for d=2!!     
   wx<-seq(0,0,length=n*length(xx[,1]))
   wx[Eff]<-quartic.v(t(t(Z)/H)[Eff,],d)
   W.x<-matrix(wx,nrow=length(xx[,1]),ncol=n,byrow=TRUE)

   
   #for (i in 1:length(xx[,1])) {        
   #     z <- x.boot - t(ll*xx[i,])
   #     wx<-seq(0,0,length=n)
   #     effect<- which(abs(t(t(z)/H))[,1]<1 & abs(t(t(z)/H))[,2]<1)    ## Only for d=2!!     
   #    wx[effect] <- quartic.v(t(t(z)/H)[effect,],d)
       # for(k in effect){    
     #    if(prod(as.numeric(z[k,]/H<=1))==1){    
       # wx[k] <- ker(z[k,]/H,d)
       # }#}
   #     W.x[i,]<-wx
   # }   
   dev<-W.x%*%psi.sample.boot
   stoch.dev[,b]<-1/(n*prod(H))*RX%*%dev # normalized deviation between theta star and theta hat   
   # print(b)
   }  #)
    
   #for(l in 1:length(xx[,1])){
   #   stoch.dev[l,]<-stoch.dev[l,]-mean(stoch.dev[l,])
   #}   
   stoch.dev<-stoch.dev-matrix(rowMeans(stoch.dev)[rep(seq_len(length(fv)),each=B)],nrow=length(xx[,1]),ncol=B,byrow=TRUE)
   max.dev.boot<-colMaxs(abs(stoch.dev))  
  # Gumbel<-sqrt(2*eta)*(max.dev.boot/K_2-d_n)
   G_alpha <- quantile(max.dev.boot, probs = 1-alpha)   
   band <- sqrt(sigma.x/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*G_alpha # 1/sqrt(n*prod(H))*sqrt(sigma.x*sqrt(K_2)/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*(d_n+G_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, sigma.x=sigma.x, Fe_x=Fe_x,eps.hat=eps.hat, max.dev=max.dev.boot,G.alpha=G_alpha)
}

lcre.boot.het<-function(x, y, bwth, d, tau = 0.5, kern="quartic", B=500, alpha=0.05,xx){
if(kern=="epanech"){
K_2 <- (3/5)^d # Product kernel Epanechnikov
SIG <- diag(d)*(3/2)*(3/5)^(d-1)  # Product kernel Epanechnikov
ker <- epanech.prod
}
if(kern=="quartic"){
K_2 <- (5/7)^d # Product kernel quartic
SIG <- diag(d)*(15/7)^(5/7)^(d-1)  # Product kernel quartic
ker <- quartic.prod
}

fv <- numeric(0)
fv.i <- numeric(0)  # for estimating the conditional varepsilon^2(x)
sigma.x <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.01)   # Undersmoothing
    ll<-matrix(1,d,n)
    W.x <- matrix(0,nrow=length(xx[,1]),ncol=n)
    for (i in 1:length(xx[,1])) {
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        W.x[i,] <- wx
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    i <- 1
    k <- 1
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv.i[i] <- r@coefficients[1]
    }
   eps.hat <- y-fv.i
   phi.i<-loss.phi(eps.hat,tau=tau)
    i <- 1
    k <- 1
     for (i in 1:length(xx[,1])){
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- lm(phi.i ~ 1, weights = wx)
        sigma.x[i] <- r$coefficients[1]
    }
   
   Fe_x.objt<-npcdist(txdat=x,tydat=eps.hat,bwmethod="normal-reference",exdat=xx,eydat=seq(0,length=length(xx[,1]),by=0))   
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=d)     
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   Fe_x<-Fe_x.objt$condist
   xx.mat<-as.matrix(xx)
   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
   kap <- abs(log(prod(H)^(1/d), base = n))
   eta <- log(vol.D)+ kap*d*log(n)
   H_2 <- vol.D*prod(H/H[1])*(2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))

   sample.boot <- matrix(0,nrow=n,ncol=B)
   dev.boot <- matrix(0,nrow=length(xx[,1]),ncol=B)
   max.dev.boot <- as.vector(seq(0,length=B,by=0))
   i <- 1
   k <- 1
   l <- 1
   for(i in 1:n){
   sample.boot[i,]<-rkde(n=B,fhat=f.star[[i]])
   }
   psi.sample.boot <- matrix(0,nrow=n,ncol=B)
   for(k in 1:B){
   psi.sample.boot[,k]<-phi(sample.boot[,k],tau)
   }
   S.FX<-diag(1/sqrt(fx*sigma.x),nrow=length(fx),ncol=length(fx))  # Diagonal matrix with 1/sqrt(sigma.x*f(x))being the diagonal line
   dev.boot<-abs(sqrt(1/(n*prod(H)))*S.FX%*%W.x%*%psi.sample.boot) # normalized deviation between theta star and theta hat
   max.dev.boot<-colMaxs(dev.boot)
   Gumbel<-sqrt(2*eta)*(max.dev.boot/K_2-d_n)
   G_alpha <- quantile(Gumbel, probs = 1-alpha)   
   band <- 1/sqrt(n*prod(H))*sqrt(sigma.x*sqrt(K_2)/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*(d_n+G_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, sigma.x=sigma.x, Fe_x=Fe_x,eps.hat=eps.hat, max.dev=max.dev.boot,G.alpha=G_alpha)
}
###################################################################################
lcre.boot.hom<-function(x, y, bwth, d, tau = 0.5, kern="quartic", B=500, alpha=0.05,xx){
if(kern=="epanech"){
K_2 <- (3/5)^d # Product kernel Epanechnikov
SIG <- diag(d)*(3/2)*(3/5)^(d-1)  # Product kernel Epanechnikov
ker <- epanech.prod
}
if(kern=="quartic"){
K_2 <- (5/7)^d # Product kernel quartic
SIG <- diag(d)*(15/7)^(5/7)^(d-1)  # Product kernel quartic
ker <- quartic.prod
}

fv <- numeric(0)
fv.i <- numeric(0)  # for estimating the conditional varepsilon^2(x)
sigma.x <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.01)   # Undersmoothing
    ll<-matrix(1,d,n)
    W.x <- matrix(0,nrow=length(xx[,1]),ncol=n)
    for (i in 1:length(xx[,1])) {
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        W.x[i,] <- wx
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    i <- 1
    k <- 1
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv.i[i] <- r@coefficients[1]
    }
   eps.hat <- y-fv.i
   phi.i<-loss.phi(eps.hat,tau=tau)
    i <- 1
    k <- 1
     for (i in 1:length(xx[,1])){
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- lm(phi.i ~ 1, weights = wx)
        sigma.x[i] <- r$coefficients[1]
    }
   
   Fe_x.objt<-npcdist(txdat=x,tydat=eps.hat,bwmethod="normal-reference",exdat=xx,eydat=seq(0,length=length(xx[,1]),by=0))
    
   i<-1
   k<-1
   
   f.star<-kde(x=eps.hat,h= Fe_x.objt$ybw,binned=TRUE)          
         
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=d)     
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   Fe_x<-Fe_x.objt$condist
   xx.mat<-as.matrix(xx)
   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
   kap <- abs(log(prod(H)^(1/d), base = n))
   eta <- log(vol.D)+ kap*d*log(n)
   H_2 <- vol.D*prod(H/H[1])*(2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))

   dev.boot <- matrix(0,nrow=length(xx[,1]),ncol=B)
   max.dev.boot <- as.vector(seq(0,length=B,by=0))
   i <- 1
   k <- 1
   l <- 1
   sample.boot<-matrix(rkde(n=n*B,fhat=f.star),nrow=n,ncol=B)
   
   psi.sample.boot <- matrix(0,nrow=n,ncol=B)
   for(k in 1:B){
   psi.sample.boot[,k]<-phi(sample.boot[,k],tau)
   }
   S.FX<-diag(1/sqrt(fx*sigma.x),nrow=length(fx),ncol=length(fx))  # Diagonal matrix with 1/sqrt(sigma.x*f(x))being the diagonal line
   dev.boot<-abs(sqrt(1/(n*prod(H)))*S.FX%*%W.x%*%psi.sample.boot) # normalized deviation between theta star and theta hat
   max.dev.boot<-colMaxs(dev.boot)
   Gumbel<-sqrt(2*eta)*(max.dev.boot/K_2-d_n)
   G_alpha <- quantile(Gumbel, probs = 1-alpha)   
   band <- 1/sqrt(n*prod(H))*sqrt(sigma.x*sqrt(K_2)/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*(d_n+G_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, sigma.x=sigma.x, Fe_x=Fe_x,eps.hat=eps.hat, max.dev=max.dev.boot,G.alpha=G_alpha)
}

######################################################################################
lcre.conf.emp<-function(x, y, bwth, d, tau = 0.5, kern="quartic", alpha=0.05,xx){
if(kern=="epanech"){
K_2 <- (3/5)^d # Product kernel Epanechnikov
SIG <- diag(d)*(3/2)*(3/5)^(d-1)  # Product kernel Epanechnikov
ker <- epanech.prod
}
if(kern=="quartic"){
K_2 <- (5/7)^d # Product kernel quartic
SIG <- diag(d)*(15/7)^(5/7)^(d-1)  # Product kernel quartic
ker <- quartic.prod
}

fv <- numeric(0)
fv.i <- numeric(0)  # for estimating the conditional varepsilon^2(x)
sigma.x <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.05)   # Undersmoothing
    for (i in 1:length(xx[,1])) {
        ll<-matrix(1,d,n)
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv[i] <- r@coefficients[1]
    }
    i <- 1
    k <- 1
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        wx[which(wx==0)]<-1.819/10^{12}
        r <- vglm(y ~ 1, weights = wx,  amlnormal(w.aml = tau/(1-tau)))
        fv.i[i] <- r@coefficients[1]
    }
   eps.hat <- y-fv.i
   phi.i<-loss.phi(eps.hat,tau=tau)
    i <- 1
    k <- 1
     for (i in 1:length(xx[,1])){
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- lm(phi.i ~ 1, weights = wx)
        sigma.x[i] <- r$coefficients[1]
    }
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=d)       
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   Fe_x.objt<-npcdist(fe_x.bw,txdat=x,tydat=eps.hat,bwmethod="normal-reference",exdat=xx,eydat=seq(0,length=length(xx[,1]),by=0))
   Fe_x<-Fe_x.objt$condist
   xx.mat<-as.matrix(xx)
   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
   kap <- abs(log(prod(H)^(1/d), base = n))
   eta <- log(vol.D)+ kap*d*log(n)
   c_alpha <- log(2)-log(abs(log(1-alpha))) 
   H_2 <- vol.D*prod(H/H[1])*(2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))
   band <- 1/sqrt(n*prod(H))*sqrt(sigma.x*sqrt(K_2)/fx)*((-2)*(Fe_x*(2*tau-1)-tau))^(-1)*(d_n+c_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, sigma.x=sigma.x, Fe_x=Fe_x)
}