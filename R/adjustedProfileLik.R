adjustedProfileLik <- function(dispersion, y, design, offset, adjust=TRUE)
# tagwise Cox-Reid adjusted profile likelihoods for the dispersion
# dispersion can be scalar or tagwise vector
# y is matrix: rows are genes/tags/transcripts, columns are samples/libraries
# offset is matrix of the same dimensions as y
# Yunshun Chen, Gordon Smyth, Aaron Lun
# Created June 2010. Last modified 21 Aug 2012.
{
	if(any(dim(y)!=dim(offset))) offset <- expandAsMatrix(offset,dim(y))
	ntags <- nrow(y)
	nlibs <- ncol(y)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

#	Fit tagwise linear models. This is actually the most time-consuming
#	operation that I can see for this function.
	fit <- glmFit(y,design=design,dispersion=dispersion,offset=offset,prior.count=0)

#	Compute log-likelihood.
	mu <- fit$fitted
	if(dispersion[1] == 0){
		loglik <- rowSums(dpois(y,lambda=mu,log = TRUE))
	} else {
		loglik <- rowSums(dnbinom(y,size=1/dispersion,mu=mu,log = TRUE))
	}
	if (!adjust) {
		return(loglik)
	}
		
#	Calculating the Cox-Reid adjustment.
	if(ncol(design)==1) {
		D <- rowSums(mu/(1+mu*dispersion))
		cr <- 0.5*log(abs(D))
	} else {
		W <- mu/(1+dispersion*mu)

#	Checking type, otherwise the C++ code will complain.
#	Note the use of a transposed matrix for easy row access in column-major format.
		if (!is.double(W)) storage.mode(W)<-"double"
		if (!is.double(design)) storage.mode(design)<-"double"
		cr <- .Call("R_cr_adjust", t(W), design, nrow(design), PACKAGE="edgeR")
		if (is.character(cr)) { stop(cr) }
	}
 
	return(loglik - cr)
}

