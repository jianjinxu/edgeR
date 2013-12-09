mglmOneGroup <- function(y,dispersion=0,offset=0,weights=NULL,maxit=50,tol=1e-10,verbose=FALSE)
#	Fit single-group negative-binomial glm
#	Aaron Lun and Gordon Smyth
#	18 Aug 2010. Last modified 22 Nov 2013.
{
#	Check y
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	if(any(y<0)) stop("y must be non-negative")
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check dispersion
	dispersion <- as.vector(dispersion)
	if(typeof(dispersion) != "double") stop("dispersion not floating point number")
	if(any(dispersion<0)) stop("dispersion must be non-negative")

#	Check offset
	if(typeof(offset) != "double") stop("offset not floating point number")

#	All-Poisson special case
	if(all(dispersion==0) && is.null(weights)) {
		N <- exp(offset)
		if(is.null(dim(N)))
			m <- mean(N)
		else
			m <- .rowMeans(N,ntags,nlibs)
	    return(log(.rowMeans(y/m, ntags, nlibs)))
	}

#	Check weights
	if(is.null(weights)) weights=1
	if(typeof(weights) == "integer") storage.mode(weights) <- "double"
	if(typeof(weights) != "double") stop("weights is non-numeric")

#	Expansions to full dimensions
	dispersion <- rep(dispersion,length=ntags)
	offset <- expandAsMatrix(offset,dim(y))
	weights <- expandAsMatrix(weights,dim(y))

#	Fisher scoring iteration.
#	Matrices are transposed so that values for each tag are in consecutive memory locations in C
	output <- .Call("R_one_group", ntags, nlibs, y, dispersion, offset, weights, maxit, tol, PACKAGE="edgeR")

#	Check error condition
	if(is.character(output)) stop(output)

#	Convergence achieved for all tags?
	if(verbose) if (any(!output[[2]])) warning(paste("max iteractions exceeded for", sum(!output[[2]]), "tags", sep=" "))

	output[[1]]
}
