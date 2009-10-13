condLogLikDerSize<-function(y,r,der=1) 
# Calculate derivatives of the conditional log-likelihood function l_g{r}
# with respect to r=1/phi (phi is the dispersion parameter)
# der is derivative (0th deriv is the function)
# For a single group of replicate libraries, all of the same total size
# Written by Mark Robinson
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
	if (der==1) {
		ll<-rowSums(digamma(y+r)) + n*digamma(n*r) - n*digamma(t+n*r) - n*digamma(r)
	} else if(der==2) {
		ll<-rowSums(trigamma(y+r)) + n^2*trigamma(n*r) - n^2*trigamma(t+n*r) - n*trigamma(r)
	} else if(der==0) {
		ll<-rowSums(lgamma(y+r)) + lgamma(n*r) - lgamma(t+n*r) - n*lgamma(r)
	}
	ll
}

