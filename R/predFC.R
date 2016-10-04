predFC <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=NULL,weights=NULL,...) 
UseMethod("predFC")

predFC.DGEList <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=NULL,weights=NULL,...)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) {
		dispersion <- 0
		message("dispersion set to zero")
	}
	predFC.default(y=y$counts,design=design,prior.count=prior.count,offset=offset,dispersion=dispersion,weights=weights,...)
}

predFC.default <- function(y,design=NULL,prior.count=0.125,offset=NULL,dispersion=0,weights=NULL,...)
#	Shrink log-fold-changes towards zero by augmenting data counts
#	Gordon Smyth and Belinda Phipson
#	17 Aug 2011.  Last modified 3 Oct 2016.
{
#	Add prior counts in proportion to library sizes
	out <- addPriorCount(y, offset=offset, prior.count=prior.count)

#	Check design
	if(is.null(design)) {
		warning("Behaviour of predFC with design=NULL is scheduled to be deprecated April 2014. Use cpm() instead.",call.=FALSE)
		return(cpm(y,lib.size=exp(offset),log=TRUE,prior.count=prior.count))
	} else
		design <- as.matrix(design)

#	Return matrix of coefficients on log2 scale
	g <- glmFit(out$y,design,offset=out$offset,dispersion=dispersion,prior.count=0,weights=weights,...)
	g$coefficients/log(2)
}

