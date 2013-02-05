mglmOneGroup <- function(y,dispersion=0,offset=0,maxit=50,tol=1e-10)
#	Fit null (single-group) negative-binomial glm with log-link to DGE data
#  Aaron Lun and Gordon Smyth
#	18 Aug 2010. Last modified 19 October 2012.
{
#	Check input values for y
	y <- as.matrix(y)
	if(any(y<0)) stop("y must be non-negative")
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check input values for dispersion
	if(any(dispersion<0)) stop("dispersion must be non-negative")


#	All-Poisson special case
	N <- exp(offset)
	if(all(dispersion==0)) {
		if(is.null(dim(N)))
			m <- mean(N)
		else
			m <- .rowMeans(N,ntags,nlibs)
		return(log(.rowMeans(y/m,ntags,nlibs)))
	}

#	Expanding the offset and dispersion values.
	dispersion <- rep(dispersion,length=ntags)
	offset <- expandAsMatrix(offset,dim(y))

#	Checking type for entry into C++ code.
	if (!is.double(dispersion)) storage.mode(dispersion)<-"double"
	if (!is.double(offset)) storage.mode(offset)<-"double"
	stopifnot(is.numeric(y));

#	Fisher scoring iteration. Matrices are transposed due to column major storage - thus, each column
#	of the transposed matrix maps to a row of the original for easy access.
	output<-.Call("R_one_group", ntags, nlibs, t(y), dispersion, t(offset), maxit, tol, PACKAGE="edgeR")
	if (is.character(output) ) { stop(output) }
	if (any(!output[[2]])) warning(paste("max iteractions exceeded for", sum(!output[[2]]), "tags", sep=" "))

	output[[1]]
}
