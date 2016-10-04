makeCompressedMatrix <- function(x, byrow=TRUE) 
# Coerces a NULL, scalar, vector or matrix to a compressed matrix,
# Determines whether the rows or columns are intended to be 
# repeated, and stores this in the attributes.
#
# written by Aaron Lun
# created 24 September 2016
# last modified 27 September 2016
{
	if (is.matrix(x)) {
		if (is(x, "compressedMatrix")) {
			return(x)
		}
		repeat.row <- repeat.col <- FALSE
	} else if (length(x)==1L) {
		repeat.row <- repeat.col <- TRUE
		x <- matrix(x)
	} else {
		if (!byrow) {
			x <- cbind(x)
			repeat.row <- FALSE
			repeat.col <- TRUE
		} else {
			x <- rbind(x)
			repeat.row <- TRUE
			repeat.col <- FALSE
		}
	}
	class(x) <- "compressedMatrix"
	attributes(x)$repeat.row <- repeat.row
	attributes(x)$repeat.col <- repeat.col
	return(x)
}

`[.compressedMatrix` <- function(x, i, j, ...)
# A wrapper function to easily subset a makeCompressedMatrix object.
#
# written by Aaron Lun
# created 24 September 2016
{
	row.status <- attributes(x)$repeat.row
	col.status <- attributes(x)$repeat.col
	oldclass <- class(x)
	x <- unclass(x)
	if (!row.status && !missing(i)) {
		x <- x[i,,drop=FALSE]
	}
	if (!col.status && !missing(j)) {
		x <- x[,j,drop=FALSE]
	}
	class(x) <- oldclass
	attributes(x)$repeat.row <- row.status
	attributes(x)$repeat.col <- col.status
	return(x)
}

as.matrix.compressedMatrix <- function(x, ...) 
# Getting rid of odds and ends.
#
# written by Aaron Lun
# created 26 September 2016
{
	attributes(x)$repeat.row <- NULL
	attributes(x)$repeat.col <- NULL
	unclass(x)
}

.addCompressedMatrices <- function(x, y) 
# A function that performs addition of two compressedMatrix objects,
# in a manner that best preserves memory usage.
#
# written by Aaron Lun
# created 26 September 2016
# last modified 27 September 2016
{
	if (!is(x, "compressedMatrix") || !is(y, "compressedMatrix")) {
		stop("only two compressedMatrix objects can be added")
	}
	dims <- pmax(dim(x), dim(y))
	out <- .Call(.cR_add_repeat_matrices, x, y, dims[1], dims[2])
	if (is.character(out)) stop(out)

	summed <- out[[1]]
	class(summed) <- class(x)
	attributes(summed)$repeat.row <- out[[2]]
	attributes(summed)$repeat.col <- out[[3]]
	return(summed)
}

.compressOffsets <- function(y, offset, lib.size=NULL) 
# A convenience function to avoid repeatedly having to write the code below.
# If provided, offsets take precedence over the library size.
# If neither are provided, library sizes are automatically computed
# as the sum of counts in the count matrix 'y'.
# A prefunctory check for finite values is performed in the C++ code.
# If 'offset' is already of the compressedMatrix class, then 
# we assume it's already gone through this once so we don't do it again.
{
	if (is(offset, "compressedMatrix")) {
		return(offset)
	}

	if (is.null(offset)) {
		if (is.null(lib.size)) lib.size <- colSums(y)
		offset <- log(lib.size)
	}
	if (!is.double(offset)) storage.mode(offset) <- "double"
	offset <- makeCompressedMatrix(offset, byrow=TRUE)

	err <- .Call(.cR_check_finite, offset, "offsets")
	if (is.character(err)) stop(err) 
	return(offset)
}

.compressWeights <- function(weights=NULL) 
# A convenience function to avoid repeatedly having to write the code below.
# All weights default to 1 if not specified.
# A prefunctory check for finite, positive values is performed in the C++ code.
# If 'weights' is already a compressedMatrix, then we assume it's 
# already gone through this and don't do it again.
{
	if (is(weights, "compressedMatrix")) {
		return(weights)
	}

	if (is.null(weights)) weights <- 1
	if (!is.double(weights)) storage.mode(weights) <- "double"
	weights <- makeCompressedMatrix(weights, byrow=TRUE)

	err <- .Call(.cR_check_positive, weights, "weights")
	if (is.character(err)) stop(err)
	return(weights)
}

.compressPrior <- function(prior.count) 
# Again for the prior counts, checking for non-negative finite values.
# Skipping the check if it's already a compressedMatrix object.
{
	if (is(prior.count, "compressedMatrix")) {
		return(prior.count)
	}
			
	if(!is.double(prior.count)) storage.mode(prior.count) <- "double"
	prior.count <- makeCompressedMatrix(prior.count, byrow=FALSE)
	err <- .Call(.cR_check_nonnegative, prior.count, "prior counts")
	if (is.character(err)) stop(err)
	return(prior.count)
}

.compressDispersions <- function(dispersion) 
# Again for the dispersions, checking for non-negative finite values.
# Skipping the check if it's already a compressedMatrix object.
{
	if (is(dispersion, "compressedMatrix")) {
		return(dispersion)
	}
			
	if(!is.double(dispersion)) storage.mode(dispersion) <- "double"
	dispersion <- makeCompressedMatrix(dispersion, byrow=FALSE)
	err <- .Call(.cR_check_nonnegative, dispersion, "dispersions")
	if (is.character(err)) stop(err)
	return(dispersion)
}

