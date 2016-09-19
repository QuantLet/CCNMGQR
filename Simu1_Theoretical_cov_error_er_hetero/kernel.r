# epanechnikov kernel
epanech.prod<-function(x,d){              # x: d*1 vector
pp<-1
for(j in 1:d){
pp<-as.numeric(abs(x[j])<1)*3/4*(1-x[j]^2)*pp
}
return(pp)
}

epanech.de.prod<-function(x,d){              # x: d*1 vector
pp<-1
for(j in 1:d){
pp<-as.numeric(abs(x[j])<1)*(-3/2)*x[j]*pp
}
return(pp)
}

############################################################################

quartic.prod<-function(x,d){              # x: d*1 vector
pp<-1
for(j in 1:d){
pp<-as.numeric(abs(x[j])<1)*15/16*(1-x[j]^2)^2*pp
}
return(pp)
}

quartic.de.prod<-function(x,d){              # x: d*1 vector
pp<-1
for(j in 1:d){
pp<-as.numeric(abs(x[j])<1)*(-15/4)*x[j]*(1-x[j]^2)*pp
}
return(pp)
}


quartic.v<-function(x,d){              # x: d*1 vector
pp<-numeric(0)
pp<-rowProds(15/16*(1-as.matrix(x)^2)^2)
return(pp)
}

 
###########################################################################3
gaussian.prod<-function(x,d){
if(d>1){
pp<-rowProds(dnorm(x,mean=0,sd=1))
}
else{
pp<-dnorm(x)
}
return(pp)
}

gaussian.de.prod<-function(x,d){
pp<-1
for(j in 1:d){
pp<-(-x[j])*dnorm(x[j],mean=0,sd=1)*pp
}
return(pp)
}