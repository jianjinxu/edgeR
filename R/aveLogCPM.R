aveLogCPM <- function(y, ...)
UseMethod("aveLogCPM")

aveLogCPM.DGEList <- function(y, normalized.lib.sizes=TRUE, prior.count=2, dispersion=NULL, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	Created 11 March 2013.  Last modified 24 November 2013.
{
#	Library sizes should be stored in y but are sometimes missing
	lib.size <- y$samples$lib.size
	if(is.null(lib.size)) lib.size <- colSums(y$counts)

#	Normalization factors should be stored in y but are sometimes missing
	if(normalized.lib.sizes) {
		nf <- y$samples$norm.factors
		if(!is.null(y$samples$norm.factors)) lib.size <- lib.size*nf
	}

#	Dispersion supplied as argument over-rules value in object
#	Should trended.dispersion or tagwise.dispersion be used instead of common.dispersion if available?
	if(is.null(dispersion)) dispersion <- y$common.dispersion

	aveLogCPM(y$counts,lib.size=lib.size,prior.count=prior.count,dispersion=dispersion,weights=y$weights)
}

aveLogCPM.DGEGLM <- function(y, prior.count=2, dispersion=NULL, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	Created 11 March 2013.  Last modified 24 Nov 2013.
{
#	Dispersion supplied as argument over-rules value in object
	if(is.null(dispersion)) dispersion <- y$dispersion

	aveLogCPM(y$counts,offset=y$offset,prior.count=prior.count,dispersion=dispersion,weights=y$weights)
}

aveLogCPM.default <- function(y,lib.size=NULL,offset=NULL,prior.count=2,dispersion=NULL,weights=NULL, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	Created 25 Aug 2012. Last modified 12 Nov 2013.
{
#	Check y
	y <- as.matrix(y)
	if(any(y<0)) stop("y must be non-negative")

#	Check prior.count
	if(prior.count<0) prior.count <- 0

#	Check dispersion
	if(is.null(dispersion)) dispersion <- 0.05
	isna <- is.na(dispersion)
	if(all(isna)) dispersion <- 0.05
	if(any(isna)) dispersion[isna] <- mean(dispersion,na.rm=TRUE)

#	Check lib.size and offset.
#	If offset is provided, it takes precedence over lib.size.
#	However it must have a similar average to log(lib.size)
#	for the results to be meaningful as logCPM values
	if(is.null(offset)) {
		if(is.null(lib.size)) lib.size <- colSums(y)
	} else {
		lib.size <- exp(offset)
	}
	mean.lib.size <- mean(lib.size)

#	Special case when all counts are zero
	if(mean.lib.size==0) {
		abundance <- rep(-log(nrow(y)),nrow(y))
		return( (abundance+log(1e6)) / log(2) )
	}

#	Scale prior counts to preserve fold changes
	prior.count.scaled <- lib.size/mean.lib.size * prior.count

#	Add double prior counts to library sizes
	offset <- log(lib.size+2*prior.count.scaled)

#	Add prior counts to y
	if(is.null(dim(prior.count.scaled))) prior.count.scaled <- matrix(1,nrow(y),1) %*% prior.count.scaled
	y <- y+prior.count.scaled

	abundance <- mglmOneGroup(y,dispersion=dispersion,offset=offset,weights=weights)
	(abundance+log(1e6)) / log(2)
}
