predFC <- function(y,design,prior.count.total=0.5,offset=NULL,dispersion=NULL) 
UseMethod("predFC")

predFC.DGEList <- function(y,design,prior.count.total=0.5,offset=NULL,dispersion=NULL)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) {
		dispersion <- 0
		message("dispersion set to zero")
	}
	predFC(y=y$counts,design=design,prior.count.total=prior.count.total,offset=offset,dispersion=dispersion)
}

predFC.default <- function(y,design,prior.count.total=0.5,offset=log(colSums(y)),dispersion=0)
#	Shrink glm estimates by augmenting data counts towards a constant
#	17 Aug 2011.  Last modified 1 May 2012.
{
	y <- as.matrix(y)
	if(missing(design)) stop("design must be set")
	ngenes <- nrow(y)
	nsamples <- ncol(y)
	pfc <- rep(0,ngenes)
	names(pfc) <- rownames(y)
	offset <- expandAsMatrix(offset,dim(y))
#	dispersion <- rep(dispersion,ngenes)

#	Add one to rowsum, in proportion to library sizes
	lib.size <- exp(offset)
	total.lib.size <- rowSums(lib.size)
	proportion <- lib.size/total.lib.size
	y.augmented <- y+proportion*prior.count.total

#	Adjust offsets to keep overall mean similar
#	total.count <- rowSums(y)
#	offset.augmented <- offset+log((total.count+1)/pmax(total.count,0.5))

   g <- glmFit(y.augmented,design,offset=offset,dispersion=dispersion,prior.count.total=0)
   g$coefficients
}

