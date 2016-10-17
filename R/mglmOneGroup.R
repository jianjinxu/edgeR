mglmOneGroup <- function(y,dispersion=0,offset=0,weights=NULL,maxit=50,tol=1e-10,verbose=FALSE,coef.start=NULL)
#	Fit single-group negative-binomial glm
#	Aaron Lun and Gordon Smyth
#	18 Aug 2010. Last modified 03 Oct 2016.
{
#	Check y
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	.isAllZero(y)

#	Check dispersion
	dispersion <- .compressDispersions(dispersion)

#	Check offset
	offset <- .compressOffsets(y, offset=offset)

#	Check starting values
	if (is.null(coef.start)) coef.start <- NA_real_
	if (!is.double(coef.start)) storage.mode(coef.start) <- "double"
	coef.start <- makeCompressedMatrix(coef.start, byrow=FALSE)

#	Check weights
	weights <- .compressWeights(weights)

#	Fisher scoring iteration.
	output <- .Call(.cR_one_group, y, dispersion, offset, weights, maxit, tol, coef.start)

#	Check error condition
	if(is.character(output)) stop(output)

#	Convergence achieved for all tags?
	if(verbose) if (any(!output[[2]])) warning(paste("max iteractions exceeded for", sum(!output[[2]]), "tags", sep=" "))

	output[[1]]
}
