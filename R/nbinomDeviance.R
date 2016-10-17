nbinomDeviance <- function(y,mean,dispersion=0,weights=NULL)
#	Negative binomial residual deviance
#	y is a matrix and a deviance is computed for each row
#	A vector y is taken to be a matrix with one row.
#	Original version 23 November 2010.
#	Last modified 03 Oct 2016.
{
#	Convert vector to row matrix
	if(!is.matrix(y)) y <- array(y, c(1L,length(y)), if(!is.null(names(y))) list(NULL,names(y)))
	out <- .compute_nbdeviance(y=y, mean=mean, dispersion=dispersion, weights=weights, dosum=TRUE)
	names(out) <- rownames(y)
	out
}

.compute_nbdeviance <- function(y, mean, dispersion, weights, dosum) {
#	Check y. May be matrix or vector.
	if(!is.matrix(y)) y <- matrix(y)
	nobs <- length(y)

#	Check mean
	if(!is.matrix(mean)) mean <- matrix(mean)
	if(!is.double(mean)) storage.mode(mean) <- "double"
	if(length(mean)<nobs) stop("mean should have same dimensions as y")

#	Check dispersion (can be tagwise (rowwise) or observation-wise).
	dispersion <- .compressDispersions(dispersion)

#	Check weights.
	weights <- .compressWeights(weights)

#	Computing unit deviance or residual deviance per gene, depending on 'dosum'.
	d <- .Call(.cR_compute_nbdev, y, mean, dispersion, weights, as.logical(dosum))
	if(is.character(d)) stop(d) 

	return(d)
}

nbinomUnitDeviance <- function(y,mean,dispersion=0) 
#	Unit deviance for the nbinom distribution.
{
	out <- .compute_nbdeviance(y=y, mean=mean, dispersion=dispersion, weights=NULL, dosum=FALSE)
	dimnames(out) <- dimnames(y)
	return(out)
}
