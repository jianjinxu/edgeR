condLogLikDerDelta<-function(y,delta,grid=TRUE,der=1,doSum=TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to calculate derivatives of the log-likelihood function, where delta is 1/(1+r)
# delta is 1/(1+r) below ... 1/delta-1=r
# der is derivative (0th deriv is the function)
# if grid=T, calculates the likelihood (derivatives) at a grid of deltas
# if grid=F, length(delta)=nrow(y) and a different delta is used for each row	
{
	if (is.vector(y)) {
		n<-length(y)
		t<-sum(y)
		g<-1    
		y<-matrix(y,nrow=1)
	} else {
		t<-rowSums(y,na.rm=TRUE)
		n<-rowSums(!is.na(y))
		g<-dim(y)[1]
	}
	if (grid==TRUE) {
		logliks<-matrix(,nrow=g,ncol=length(delta))
		for (i in seq(along=delta)) {
			d<-delta[i]
			r<-(1/d)-1
			if (der==1) {
				ll<-condLogLikDerSize(y,r,der=1)*(-d^(-2))
			} else if(der==2) {
				ll<-condLogLikDerSize(y,r,der=1)*2*(d^(-3))+condLogLikDerSize(y,r,der=2)*(d^(-4))
			} else if(der==0) {
				ll<-condLogLikDerSize(y,r,der=0)
			}
			logliks[,i]<-ll
		}
		if (doSum) {
			return(colSums(logliks))
		} else {
			return(logliks)
		}
	} else {
		if( !(length(delta)==1 | length(delta)==nrow(y)) ) {
			stop("When grid=FALSE, delta must be length 1 or nrow(y)")
		}
		r<-(1/delta)-1
		if (der==1) {
			ll<-condLogLikDerSize(y,r,der=1)*(-delta^(-2))
		} else if(der==2) {
			ll<-condLogLikDerSize(y,r,der=1)*2*(delta^(-3))+condLogLikDerSize(y,r,der=2)*(delta^(-4))
		} else if(der==0) {
			ll<-condLogLikDerSize(y,r,der=0)
		}
		if (doSum)
			return(sum(ll))
		else 
			return(ll)
	}
}

