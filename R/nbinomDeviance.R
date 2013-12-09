nbinomDeviance <- function(y,mean,dispersion=0,weights=NULL)
#	Negative binomial residual deviance
#	y is a matrix, and a deviance is computed to be a row
#	A vector y is taken to be a matrix with one row.
#	Original version 23 November 2010.
#	Last modified 9 Dec 2013.
{
#	Convert vector to row matrix
	if(!is.matrix(y)) y <- array(y, c(1L,length(y)), if(!is.null(names(y))) list(NULL,names(y)))

	d <- nbinomUnitDeviance(y=y,mean=mean,dispersion=dispersion)
	if(!is.null(weights)) d <- weights*d
	rowSums(d)
}


nbinomUnitDeviance <- function(y,mean,dispersion=0) 
#	Unit deviance for the nbinom distribution.
{
#	Check y. May be matrix or vector.
	if (!is.double(y)) storage.mode(y) <- "double"
	ntags <- NROW(y)
	nobs <- length(y)

#	Check mean
	if (!is.double(mean)) storage.mode(mean) <- "double"
	if(length(mean)<nobs) stop("mean should have same dimensions as y")

#	Check dispersion.
#	Can be tagwise (rowwise) or observation-wise.
	if (!is.double(dispersion)) dispersion <- "double"
	lend <- length(dispersion)
	if(lend < ntags) dispersion <- rep_len(dispersion, length.out=ntags)
	if(lend > ntags && lend < nobs) dispersion <- rep_len(dispersion, length.out=nobs)

	out <- .Call("R_compute_nbdev", y=y, mu=mean, phi=dispersion, PACKAGE="edgeR")

#	Check error status
	if (is.character(out)) stop(out)

	y[] <- out
	return(y)
}
