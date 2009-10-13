logLikDerP<-function(p,y,lib.size,r,der=0) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function calculating derivatives of the log-likelihood w.r.t. p (assuming r known)
{
	if (is.vector(y)) {
		y<-matrix(y,nrow=1)
	}
	n<-outer(rep(1,nrow(y)),lib.size)
	pm<-outer(p,rep(1,ncol(y)))
	if (length(r)==nrow(y)) {
		rmm<-outer(r,rep(1,ncol(y)))
	} else {
		rmm<-matrix(r,nrow=nrow(y),ncol=ncol(y))
	}
	if (der==1) {
		ret<-(1/p)*rowSums(rmm*(y-n*pm)/(n*pm+rmm),na.rm=TRUE)
		ret[p==0]<-0
	} else if (der==0) {
		ret<-rowSums(-(rmm+y)*log(rmm+n*pm)+y*log(n*pm),na.rm=TRUE) 
		ret[p==0]<-1
	} else if (der==2) {
		ret<-(1/p)*rowSums( -(rmm/p)*(y-n*pm)/(n*pm+rmm) - rmm*n*(y+rmm)/((n*pm+rmm)^2),na.rm=TRUE)
		ret[p==0]<-1
	}
	return(ret)
}

