

.interpolateHelperNew<-function(mu,p,r,d,N,lib.size,verbose=TRUE,every=1000) {
   if(length(r)==1) { r<-matrix(r,nrow=nrow(mu),ncol=ncol(mu)) }
   if(length(r)==nrow(mu)) { r<-outer(r,rep(1,ncol(mu))) }
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   psudo<-matrix(,nrow=nrow(p),ncol=ncol(p))
   for (i in 1:nrow(psudo)) {
     mx<-max(qnbinom(p[i,],size=r[i,],mu=mu[i,]))+1
     if(is.infinite(mx))
       mx<-max(d[i,])*1.5
     for (j in 1:ncol(psudo)) {
       if(lib.size[j] < N) {     
         # the count will increase ------------------------------
         for(k in d[i,j]:mx) { 
           left<-pnbinom(k-1,size=2,mu=mu[i,j])+dnbinom(k,size=2,mu=mu[i,j])/2 
           right<-pnbinom(k,size=2,mu=mu[i,j])+dnbinom(k+1,size=2,mu=mu[i,j])/2 
           if( right > p[i,j] ) {
             psudo[i,j]<-approx(c(left,right),c(d[i,j],d[i,j]+1),xout=p[i,j])$y
             break
           }
         }
         # ------------------------------------------------------
         #a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
         #psudo[i,j]<-approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y
       } else if (lib.size[j] > N) { # count will decrease
         # the count will increase ------------------------------
         for(k in d[i,j]:0) {  
           right<-pnbinom(k-1,size=2,mu=mu[i,j])+dnbinom(k,size=2,mu=mu[i,j])/2
           left<-pnbinom(k-2,size=2,mu=mu[i,j])+dnbinom(k-1,size=2,mu=mu[i,j])/2 
           if(left < p[i,j] ) {
             psudo[i,j]<-approx(c(left,right),c(d[i,j]-1,d[i,j]),xout=p[i,j])$y
             break
           }
         }
         #a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
         #psudo[i,j]<-approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y
         # ------------------------------------------------------
       } else {                  # count will stay same
         psudo[i,j]<-d[i,j]
       }
     }
     if (verbose) { if(i%%every==0) cat(".") }
   }
   v<-which(is.na(psudo),arr.ind=TRUE)
   mx<-max(qnbinom(p,size=r,mu=mu))+1
   if(is.infinite(mx))
     mx<-max(d)*1.5
   if( nrow(v) > 0 ) {
     for(i in 1:nrow(v)) {
       a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
       psudo[v[i,1],v[i,2]]<-approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y
       if (verbose) cat("-")
     }
   }
   if (verbose) cat("\n")
   return(psudo)
}

