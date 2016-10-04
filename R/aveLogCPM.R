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
#	Created 25 Aug 2012. Last modified 25 Sep 2014.
{
#	Special case when all counts and library sizes are zero
	y <- as.matrix(y)
	if(.isAllZero(y)) {
		if ((is.null(lib.size) || max(lib.size)==0) && (is.null(offset) || max(offset)==-Inf)) {
			abundance <- rep(-log(nrow(y)),nrow(y))
			return( (abundance+log(1e6)) / log(2) )
		}
	}

#	Check dispersion
	if(is.null(dispersion)) dispersion <- 0.05
	isna <- is.na(dispersion) # ???
	if(all(isna)) dispersion <- 0.05
	if(any(isna)) dispersion[isna] <- mean(dispersion,na.rm=TRUE)

	dispersion <- .compressDispersions(dispersion)

#   Check weights
	weights <- .compressWeights(weights)

#   Check offsets
	offset <- .compressOffsets(y, lib.size=lib.size, offset=offset)

#   Check prior counts
	prior.count <- .compressPrior(prior.count)

#   Retrieve GLM fitting parameters
	maxit <- formals(mglmOneGroup)$maxit
	tol <- formals(mglmOneGroup)$tol

#   Calling the C++ code
	ab <- .Call(.cR_ave_log_cpm, y, offset, prior.count, dispersion, weights, maxit, tol)
	if (is.character(ab)) stop(ab)

	return(ab)
}

