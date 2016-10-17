designAsFactor <- function(design)
#	Construct a factor from the unique rows of a matrix
#	Gordon Smyth
#	11 March 2011.  Last modified 19 March 2011.
{
	design <- as.matrix(design)
	z <- (exp(1)+pi)/5
	g <- factor(rowMeans(design*z^(col(design)-1)))
	levels(g) <- 1:length(levels(g))
	g
}

mglmOneWay <- function(y,design=NULL,dispersion=0,offset=0,weights=NULL,maxit=50,tol=1e-10,coef.start=NULL)
#	Fit multiple negative binomial glms with log link
#	by Fisher scoring with
#	only a single explanatory factor in the model
#	Gordon Smyth
#	11 March 2011.  Last modified 03 Oct 2016.
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)
	if(is.null(design)) {
		design <- matrix(1,nlibs,1)
		group <- factor(design)
	} else {
		design <- as.matrix(design)
		group <- designAsFactor(design)
	}
	ngroups <- length(levels(group))
	stopifnot(ncol(design)==ngroups)
	groupbeta <- matrix(0,ntags,ngroups)

	offset <- .compressOffsets(y, offset=offset)
	dispersion <- .compressDispersions(dispersion)
	weights <- .compressWeights(weights)
	
	firstjofgroup <- rep(0,ngroups)
	new.start <- NULL
	for (g in seq_len(ngroups)) {
		j <- which(group==(levels(group)[g]))
		firstjofgroup[g] <- j[1]
		if (!is.null(coef.start)) { new.start <- coef.start %*% design[firstjofgroup[g],] }
		groupbeta[,g] <- mglmOneGroup(y[,j,drop=FALSE], dispersion=dispersion[,j,drop=FALSE], offset=offset[,j,drop=FALSE],
					weights=weights[,j,drop=FALSE], maxit=maxit,tol=tol, coef.start=new.start)
	}

	# Computing the fitted values from the group-wise beta's.
	mu <- .Call(.cR_get_one_way_fitted, groupbeta, offset, as.integer(group)-1L)
	if (is.character(mu)) stop(mu)

	# Reformatting the beta's to reflect the original design.
	designunique <- design[firstjofgroup,,drop=FALSE]
	groupbeta <- pmax(groupbeta,-1e8)
	beta <- t(solve(designunique,t(groupbeta)))

	list(coefficients=beta,fitted.values=mu)
}
