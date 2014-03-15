#  sumTechReps.R

sumTechReps <- function(x,ID=colnames(x),...) UseMethod("sumTechReps")
#  24 Sept 2010

sumTechReps.default <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for matrices
#	Yifang Hu and Gordon Smyth
#	Created 14 March 2014
{
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	if(is.null(ID)) stop("No sample IDs")
	t(rowsum(t(x),group=ID,reorder=FALSE,na.rm=FALSE))
}

sumTechReps.DGEList <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for matrices
#	Yifang Hu and Gordon Smyth
#	Created 14 March 2014
{
	x$common.dispersion <- x$trended.dispersion <- x$tagwise.dispersion <- NULL
	x$weights <- NULL
	d <- duplicated(ID)
	if(!any(d)) return(x)
	y <- x[,!d]
	y$counts <- sumTechReps(x$counts,ID=ID,...)
	y$samples$lib.size <- rowsum(x$samples$lib.size,group=ID,reorder=FALSE,na.rm=FALSE)
	y$samples$norm.factors <- tapply(x$samples$norm.factors,ID,mean)
	rownames(y$samples) <- colnames(y$counts)
	y
}
