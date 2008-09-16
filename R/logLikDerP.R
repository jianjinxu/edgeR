

logLikDerP<-function(p,y,lib.size,r,der=0) {
  # log-likelihood (assuming r known), derivatives w.r.t. p
  if (is.vector(y)) {
    y<-matrix(y,nrow=1)
  }
  n<-outer(rep(1,nrow(y)),lib.size)
  pm<-outer(p,rep(1,ncol(y)))
  if (length(r)==nrow(y)) {
    rm<-outer(r,rep(1,ncol(y)))
  } else {
    rm<-matrix(r,nrow=nrow(y),ncol=ncol(y))
  }
  if (der==1) {
     ret<-(1/p)*rowSums(rm*(y-n*pm)/(n*pm+rm),na.rm=TRUE)
     ret[p==0]<-0
  } else if (der==0) {
     ret<-rowSums(-(rm+y)*log(rm+n*pm)+y*log(n*pm),na.rm=TRUE) 
     ret[p==0]<-1
  } else if (der==2) {
     ret<-(1/p)*rowSums( -(rm/p)*(y-n*pm)/(n*pm+rm) - rm*n*(y+rm)/((n*pm+rm)^2),na.rm=TRUE)
     ret[p==0]<-1
  }
  return(ret)
}

