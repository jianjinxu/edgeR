aveLogCPM <- function(y, ...)
UseMethod("aveLogCPM")

aveLogCPM.DGEList <- function(y, normalized.lib.sizes=TRUE, prior.count=2, dispersion=0.05, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	11 March 2013.
{
	lib.size <- y$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors
	aveLogCPM(y$counts,lib.size=lib.size,prior.count=prior.count,dispersion=dispersion)
}

aveLogCPM.DGEGLM <- function(y, prior.count=2, dispersion=0.05, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	11 March 2013.
{
	offset <- y$offset
	if(is.matrix(offset)) offset <- colMeans(offset)
	lib.size <- exp(offset)
	aveLogCPM(y$counts,lib.size=lib.size,prior.count=prior.count,dispersion=dispersion)
}

aveLogCPM.default <- function(y,lib.size=NULL,prior.count=2,dispersion=0.05, ...)
#	log2(AveCPM)
#	Gordon Smyth
#	25 Aug 2012. Last modified 30 Sept 2013.
{
#	Check input
	y <- as.matrix(y)
	if(any(y<0)) stop("y must be non-negative")
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(prior.count<0) prior.count <- 0

	mean.lib.size <- mean(lib.size)
	if(mean.lib.size==0) {
		abundance <- rep(-log(nrow(y)),nrow(y))
	} else {
		prior.count.scaled <- lib.size/mean.lib.size * prior.count
		offset <- log(lib.size+2*prior.count.scaled)
		abundance <- mglmOneGroup(t(t(y)+prior.count.scaled),dispersion=dispersion,offset=offset)
	}
	(abundance+log(1e6)) / log(2)
}
