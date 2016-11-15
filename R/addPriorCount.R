addPriorCount <- function(y, lib.size=NULL, offset=NULL, prior.count=1) 
# Add library size-adjusted prior counts to values of 'y'.
# Also add twice the adjusted prior to the library sizes, 
# which are provided as log-transformed values in 'offset'.
#
# written by Aaron Lun
# created 26 September 2016
{
#	Check y
	if (!is.numeric(y)) stop('count matrix must be numeric')
	y <- as.matrix(y)

#	Check prior.count
	prior.count <- .compressPrior(prior.count)

#	Check lib.size and offset.
#	If offsets are provided, they must have a similar average to log(lib.size)
#	for the results to be meaningful as logCPM values
	offset <- .compressOffsets(y, lib.size=lib.size, offset=offset)

#   Adding the prior count.
	out <- .Call(.cR_add_prior_count, y, offset, prior.count)
	if (is.character(out)) stop(out)
    
	names(out) <- c("y", "offset")
	out$offset <- makeCompressedMatrix(out$offset, byrow=TRUE)
	return(out)
}

