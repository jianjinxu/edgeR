cpm <- function(x, ...)
UseMethod("cpm")

cpm.DGEList <- function(x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#	Counts per million for a DGEList
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 6 July 2015
{
	lib.size <- x$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*x$samples$norm.factors
	cpm.default(x$counts,lib.size=lib.size,log=log,prior.count=prior.count)
}

cpm.default <- function(x, lib.size=NULL, log=FALSE, prior.count=0.25, ...)
#	Counts per million for a matrix
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 03 October 2016
{
	x <- as.matrix(x)
	if (any(dim(x)==0L)) {
		return(x)
	}

	if(is.null(lib.size)) lib.size <- colSums(x)
	if(!is.double(lib.size)) storage.mode(lib.size) <- "double"
	lib.size <- makeCompressedMatrix(lib.size, byrow=TRUE)
	err <- .Call(.cR_check_positive, lib.size, "library sizes")
	if (is.character(err)) stop(err)

	# Calculating in C++ for max efficiency
	if(log) {
		prior.count <- .compressPrior(prior.count)
		out <- .Call(.cR_calculate_cpm_log, x, lib.size, prior.count)
	} else {
		out <- .Call(.cR_calculate_cpm_raw, x, lib.size)
	}
	if (is.character(out)) stop(out)

	# Cleaning up
	dimnames(out) <- dimnames(x)
	out
}
