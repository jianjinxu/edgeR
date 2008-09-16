

interpolateHelper<-function(mu,p,r,d,verbose=TRUE) {
   if(length(r)==1) { r<-matrix(r,nrow=nrow(mu),ncol=ncol(mu)) }
   if(length(r)==nrow(mu)) { r<-outer(r,rep(1,ncol(mu))) }
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   psudo<-matrix(,nrow=nrow(p),ncol=ncol(p))
   for (i in 1:nrow(psudo)) {
     mx<-max(qnbinom(p[i,],size=r[i,],mu=mu[i,]))+1
	 if(is.infinite(mx))
	   mx<-max(d[i,])*1.5
     for (j in 1:ncol(psudo)) {
       a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
       psudo[i,j]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
     }
     if (verbose) { if(i%%1000==0) cat(".") }
   }
   v<-which(is.na(psudo),arr.ind=TRUE)
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   if(is.infinite(mx))
     mx<-max(d)*1.5
   if( nrow(v) > 0 ) {
     for(i in 1:nrow(v)) {
       a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
       psudo[v[i,1],v[i,2]]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
       if (verbose) cat("-")
     }
   }
   if (verbose) cat("\n")
   return(psudo)
}

