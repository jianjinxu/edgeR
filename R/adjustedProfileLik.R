adjustedProfileLik <- function(dispersion, y, design, offset, weights=NULL, adjust=TRUE)
# tagwise Cox-Reid adjusted profile likelihoods for the dispersion
# dispersion can be scalar or tagwise vector
# y is matrix: rows are genes/tags/transcripts, columns are samples/libraries
# offset is matrix of the same dimensions as y
# Yunshun Chen, Gordon Smyth, Aaron Lun
# Created June 2010. Last modified 17 Feb 2014.
{
	if(any(dim(y)!=dim(offset))) offset <- expandAsMatrix(offset,dim(y))
	ntags <- nrow(y)
	nlibs <- ncol(y)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)
       
#	Fit tagwise linear models. This is actually the most time-consuming
#	operation that I can see for this function.
	fit <- glmFit(y,design=design,dispersion=dispersion,offset=offset,prior.count=0,weights=weights)

#	Compute log-likelihood.
	mu <- fit$fitted
	if(is.null(weights)) weights <- 1
	if(dispersion[1] == 0){
#		loglik <- rowSums(weights*dpois(y,lambda=mu,log = TRUE))
		ll <- y*log(mu)-mu-lgamma(y+1)
		ll[mu==0] <- 0
		loglik <- rowSums(weights*ll)
		
	} else {
#		loglik <- rowSums(weights*dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
		r <- 1/dispersion
		ll <- y*log(mu)-y*log(mu+r)+r*log(r)-r*log(mu+r)+lgamma(y+r)-lgamma(y+1)-lgamma(r)
		ll[mu==0] <- 0
		loglik <- rowSums(weights*ll)		
	}
	if (!adjust) {
		return(loglik)
	}
		
#	Calculating the Cox-Reid adjustment.
	if(ncol(design)==1) {
		D <- rowSums(weights*mu/(1+mu*dispersion))
		cr <- 0.5*log(abs(D))
	} else {
		W <- weights*mu/(1+dispersion*mu)

#	Checking type, otherwise the C++ code will complain.
#	Note the use of a transposed matrix for easy row access in column-major format.
		if (!is.double(W)) storage.mode(W)<-"double"
		if (!is.double(design)) storage.mode(design)<-"double"
		cr <- .Call("R_cr_adjust", t(W), design, nrow(design), PACKAGE="edgeR")
		if (is.character(cr)) { stop(cr) }
	}
 
	return(loglik - cr)
}

