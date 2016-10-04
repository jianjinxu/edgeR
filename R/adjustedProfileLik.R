adjustedProfileLik <- function(dispersion, y, design, offset, weights=NULL, adjust=TRUE, start=NULL, get.coef=FALSE)
# tagwise Cox-Reid adjusted profile likelihoods for the dispersion
# dispersion can be scalar or tagwise vector
# y is matrix: rows are genes/tags/transcripts, columns are samples/libraries
# offset is matrix of the same dimensions as y
# Yunshun Chen, Gordon Smyth, Aaron Lun
# Created June 2010. Last modified 27 July 2015.
{
#   Checking counts
	if (!is.numeric(y)) stop("counts must be numeric")
	y <- as.matrix(y)

#   Checking offsets
	if (!is.double(offset)) storage.mode(offset) <- "double"
	offset <- makeCompressedMatrix(offset, byrow=TRUE)

#   Checking dispersion
	if (!is.double(dispersion)) storage.mode(dispersion) <- "double"
	dispersion <- makeCompressedMatrix(dispersion, byrow=FALSE)

#   Checking weights
	if(is.null(weights)) weights <- 1
	if (!is.double(weights)) storage.mode(weights) <- "double"
	weights <- makeCompressedMatrix(weights, byrow=TRUE)
	  
#	Fit tagwise linear models. This is actually the most time-consuming
#	operation that I can see for this function.
	fit <- glmFit(y,design=design,dispersion=dispersion,offset=offset,prior.count=0,weights=weights,start=start)
	mu <- fit$fitted

#   Check other inputs to C++ code
	adjust <- as.logical(adjust)
	if (!is.double(design)) storage.mode(design)<-"double"

#   Compute adjusted log-likelihood
	apl <- .Call(.cR_compute_apl, y, mu, dispersion, weights, adjust, design)
	if (is.character(apl)) stop(apl)

#	Deciding what to return.
	if (get.coef) { 
		return(list(apl=apl, beta=fit$coefficients))
	} else {
		return(apl)
	}
}

