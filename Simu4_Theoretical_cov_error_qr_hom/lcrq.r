# Depend on packages: quantreg, np, matrixStats, ks
# Notice: we apply special case for kernel choice: PRODUCT kernel

lcrq<-function(x, y, H, d, tau = 0.5, xx)  # H: sequence of bandwidth, xx: a matrix m*d, x: matrix n*d, y: n*1 vector
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
        r <- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE)
        fv[i] <- r$coef[1]
    }
    list(xx = xx, fv = fv)
}

llrq<-function(x, y, bwth, d, tau = 0.5, kern="quartic", xx)  # H: sequence of bandwidth, xx: a matrix m*d, x: matrix n*d, y: n*1 vector
{
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
    dv <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw
    for (i in 1:length(xx[,1])) {
        ll<-matrix(1,d,n)
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
        fv[i] <- r$coef[1]
        dv[i] <- r$coef[2]
    }
    list(xx = xx, fv = fv,dv=dv)
}

###########################################################################
loss.psi<-function(x,tau){
phii<-numeric(0)
phii<-as.numeric(x<0)-tau
return(phii)
}
#############################################################################
lcrq.conf<-function(x, y, bwth, d, tau = 0.5, kern="quartic", alpha=0.05,xx){
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
fv.i <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw*(4*(tau*(1-tau)/(dnorm(qnorm(p=tau))^2)))^(1/5)/n^(0.05) # Yu and Jones 1998 rule of thumb (simuation with normal) + Undersmoothing
    for (i in 1:length(xx[,1])) {
        ll<-matrix(1,d,n)
        z <- x - t(ll*xx[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H,d)
        }
        r <- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE)
        fv[i] <- r$coef[1]
    }
    i <- 1
    k <- 1
    H.e <- H <- bwth$xbw*(4*(tau*(1-tau)/(dnorm(qnorm(p=tau))^2)))^(1/5) 
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H.e,d)
        }
        r <- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE)
        fv.i[i] <- r$coef[1]
    }
    eps.hat <- y-fv.i
    
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=2)       
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
   fe_x.bw<-npcdensbw(xdat=x,ydat=eps.hat,bwmethod="normal-reference",ckerorder=2)             
   fe_x.objt<-npcdens(fe_x.bw,txdat=x,tydat=eps.hat,exdat=xx,eydat=seq(0,length=length(xx[,1]),by=0))
   fy_x.bw<-fe_x.bw
   fy_x.bw$xbw<-fy_x.bw$xbw/3
   fym_x.objt<-npcdens(fy_x.bw,txdat=x,tydat=y,exdat=xx,eydat=fv)
   fe_x<-fe_x.objt$condens
   fym_x<-fym_x.objt$condens
   
   xx.mat<-as.matrix(xx)
   vol.D <- prod(colMaxs(xx.mat)-colMins(xx.mat))
   kap <- abs(log(prod(H)^(1/d), base = n))
   eta <- log(vol.D)+ kap*d*log(n)
   c_alpha <- log(2)-log(abs(log(1-alpha))) 
   H_2 <- (2*pi*K_2)^(-d/2)*sqrt(det(SIG))  
   d_n <- sqrt(2*eta)+sqrt(2*eta)^(-1)*(0.5*(d-1)*log(log(n^kap))+log(sqrt(2*pi)^(-1)*H_2*(2*d)^((d-1)/2)))
   band <- 1/sqrt(n*prod(H))*sqrt(tau*(1-tau)*sqrt(K_2)/fx)*(1/fe_x)*(d_n+c_alpha*sqrt(2*eta)^(-1))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, fe_x=fe_x, fym_x=fym_x,eps.hat=eps.hat)
}
############################################################################################################################




lcrq.boot<-function(x, y, bwth, d, tau = 0.5, kern="quartic", B=500,alpha=0.05,xx){
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
fv.i <- numeric(0)
    # dv <- xx
    n <- length(y)
    H <- bwth$xbw/n^(0.05)*(4*(tau*(1-tau)/(dnorm(qnorm(p=tau))^2)))^(1/5) # Yu and Jones 1998 rule of thumb (simuation with normal) + Undersmoothing
    ll<-matrix(1,d,n)
    W.xx.n<-matrix(0,nrow=length(xx[,1]),ncol=n)
    for (i in 1:length(xx[,1])) {        
        z <- x - t(ll*xx[i,])
        wx<-seq(0,0,length=n)
        for(k in 1:n){
         if(prod(as.numeric(z[k,]/H<=1))==1){
        wx[k] <- ker(z[k,]/H,d)
        }}
        W.xx.n[i,] <- wx   
        r <- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE)
        fv[i] <- r$coef[1]
    }   
    i <- 1
    k <- 1
    H.e <- H <- bwth$xbw*(4*(tau*(1-tau)/(dnorm(qnorm(p=tau))^2)))^(1/5)
    for (i in 1:n) {
        z <- x - t(ll*x[i,])
        wx<-numeric(0)
        for(k in 1:n){
        wx[k] <- ker(z[k,]/H.e,d)
        }
        r <- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE)
        fv.i[i] <- r$coef[1]
    }
    eps.hat <- y-fv.i
   fx.bw<-npudensbw(dat=x,bwmethod="normal-reference",ckerorder=2)       
   fx.objt<-npudens(fx.bw,tdat=x,edat=xx)
   fx<-fx.objt$dens
     
   fe_x.bw<-npcdensbw(xdat=x,ydat=eps.hat,bwmethod="normal-reference",ckerorder=2)          

   fe_x.objt<-npcdens(fe_x.bw,txdat=x,tydat=y,exdat=xx,eydat=fv)
   fe_x<-fe_x.objt$condens

   fex.bw$bw[1]<-fe_x$bw[1]*1.5
   fyx.objt<-npudens(fex.bw,tdat=cbind(y,x),edat=cbind(fv,xx))
   fyx<-fyx.objt$dens
 ######   Bootstrap loop         #####################################  
   FX<-diag(1/sqrt(fx)*(fyx/fex),nrow=length(fx),ncol=length(fx))  # Diagonal matrix with 1/sqrt(f(x))being the diagonal line

   BOOT <- cbind(x,eps.hat)  # epsilon is at the last column of this matrix
   sample.boot <- numeric(0)
   HH<-fe_x.bw$xbw   
   HH[d+1]<-fe_x.bw$ybw
   stoch.dev <- matrix(0,nrow=length(xx[,1]),ncol=B)
   max.dev.boot <- as.vector(seq(0,length=B,by=0))
   W.x.b <- matrix(0,nrow=length(xx[,1]),ncol=n)
   i <- 1
   k <- 1
   l <- 1
   for(b in 1:B){
   bb<-sample(c(1:n),size=n,replace=TRUE)   
   er<-rnorm(n*(d+1))
   E<-matrix(er,nrow=n,ncol=d+1)
   boot.data<-as.matrix(BOOT[bb,])+ E%*%diag(HH,nrow=d+1,ncol=d+1)

   x.boot <- boot.data[,1:d]
   sample.boot <- boot.data[,d+1]
   psi.sample.boot <- as.matrix(loss.psi(sample.boot,tau))
   W.x<-matrix(0,nrow=length(xx[,1]),ncol=n)
   X.boot<-matrix(t(x.boot),nrow=n*length(xx[,1]),ncol=d,byrow=TRUE)
   XX<-xx[rep(seq_len(nrow(xx)),each=n),]   
   Z<-X.boot-XX
   Eff<- which(abs(t(t(Z)/H))[,1]<1 & abs(t(t(Z)/H))[,2]<1)    ## Only for d=2!!     
   wx<-seq(0,0,length=n*length(xx[,1]))
   wx[Eff]<-quartic.v(t(t(Z)/H)[Eff,],d)
   W.x<-matrix(wx,nrow=length(xx[,1]),ncol=n,byrow=TRUE)
   dev<-W.x%*%psi.sample.boot
   stoch.dev[,b]<-(n*prod(H))^(-1)*sqrt(1/(tau*(1-tau)))*FX%*%dev
   }  #)
       
   stoch.dev<-stoch.dev-matrix(rowMeans(stoch.dev)[rep(seq_len(length(fv)),each=B)],nrow=length(xx[,1]),ncol=B,byrow=TRUE)
   max.dev.boot<-colMaxs(abs(stoch.dev))
   
   G_alpha <- quantile(max.dev.boot, probs = 1-alpha)   
   band <- sqrt(tau*(1-tau)/fx)*(1/fe_x)*G_alpha# 1/sqrt(n*prod(H))*sqrt(tau*(1-tau)*sqrt(K_2)/fx)*(1/fe_x)*(d_n+G_alpha/sqrt(2*eta))
   list(xx = xx, fv = fv,hband = fv + band, lband = fv - band, xdensbw=1/fx.bw$bw,bw = H,fx=fx, fe_x=fe_x,eps.hat=eps.hat,max.dev=max.dev.boot,G.alpha=G_alpha,boot.bw=HH)
}
