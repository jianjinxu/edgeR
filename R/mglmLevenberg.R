mglmLevenberg <- function(y, design, dispersion=0, offset=0, weights=NULL, coef.start=NULL, start.method="null", maxit=200, tol=1e-06)
#	Fit genewise negative binomial glms with log-link
#	using Levenberg damping to ensure convergence

#	R version by Gordon Smyth and Yunshun Chen
#	C++ version by Aaron Lun
#	Created 3 March 2011.  Last modified 03 Oct 2016
{
#	Check arguments
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	if(nlibs==0 || ngenes==0) stop("no data")

#	Checks for negative, NA or non-finite values in the count matrix.
	.isAllZero(y)

#	Checking the design matrix
	design <- as.matrix(design)
	if (!is.double(design)) storage.mode(design) <- "double"
	if (!all(is.finite(design))) stop("all entries of design matrix must be finite and non-missing")

#	Checking dispersions, offsets and weights
	offset <- .compressOffsets(y, offset=offset)
    dispersion <- .compressDispersions(dispersion)
	weights <- .compressWeights(weights)

#	Initializing values for the coefficients at reasonable best guess with linear models.
	if(is.null(coef.start)) {
		start.method <- match.arg(start.method, c("null","y"))
		beta <- .Call(.cR_get_levenberg_start, y, offset, dispersion, weights, design, start.method=="null")
		if (is.character(beta)) stop(beta) 
	} else {
		beta <- as.matrix(coef.start)
	}

# 	Checking arguments and calling the C++ method.
	if (!is.double(beta)) storage.mode(beta) <- "double"
	output <- .Call(.cR_levenberg, design, y, dispersion, offset, weights, beta, tol, maxit)

#	Check for error condition
	if (is.character(output)) { stop(output) }

#	Naming the output and returning it.  
	names(output) <- c("coefficients", "fitted.values", "deviance", "iter", "failed")
	colnames(output$coefficients) <- colnames(design)
	rownames(output$coefficients) <- rownames(y)
	dimnames(output$fitted.values) <- dimnames(y)
	output
}
