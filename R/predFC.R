predFC <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=NULL) 
UseMethod("predFC")

predFC.DGEList <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=NULL)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) {
		dispersion <- 0
		message("dispersion set to zero")
	}
	predFC.default(y=y$counts,design=design,prior.count=prior.count,offset=offset,dispersion=dispersion)
}

predFC.default <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=0)
#	Shrink log-fold-changes towards zero by augmenting data counts
#	Gordon Smyth and Belinda Phipson
#	17 Aug 2011.  Last modified 4 Nov 2012.
{
#	Check y
	y <- as.matrix(y)
	ngenes <- nrow(y)
	nsamples <- ncol(y)

#	Check prior.count
	if(prior.count<0) stop("prior.count should be non-negative")

#	Check offset
	if(is.null(offset)) {
		lib.size <- colSums(y)
		offset <- log(lib.size)
	} else
		lib.size <- exp(offset)

#	Check design
	if(is.null(design)) {
		warning("Behaviour of predFC with design=NULL is scheduled to be deprecated April 2014. Use cpm() instead.",call.=FALSE)
		return(cpm(y,lib.size=lib.size,log=TRUE,prior.count=prior.count))
	} else
		design <- as.matrix(design)

#	Add prior counts in proportion to library sizes
	if(is.null(dim(lib.size)))
		ave.lib.size <- mean(lib.size)
	else
		ave.lib.size <- rowMeans(lib.size)
	prior.count <- prior.count * lib.size/ave.lib.size
	lib.size <- lib.size+2*prior.count
	if(is.null(dim(prior.count))) prior.count <- matrix(prior.count,ngenes,nsamples,byrow=TRUE)
	y <- y+prior.count

#	Return matrix of coefficients on log2 scale
	g <- glmFit(y,design,offset=log(lib.size),dispersion=dispersion,prior.count=0)
	g$coefficients/log(2)
}

