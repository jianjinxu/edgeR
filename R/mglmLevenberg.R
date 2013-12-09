mglmLevenberg <- function(y, design, dispersion=0, offset=0, weights=NULL, coef.start=NULL, start.method="null", maxit=200, tol=1e-06)
#	Fit genewise negative binomial glms with log-link
#	using Levenberg damping to ensure convergence

#	R version by Gordon Smyth and Yunshun Chen
#	C++ version by Aaron Lun
#	Created 3 March 2011.  Last modified 11 July 2012
{
#	Check arguments
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	if(any(y<0)) stop("y must be non-negative")
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	if(nlibs==0 || ngenes==0) stop("no data")

	if(!( all(is.finite(y)) || all(is.finite(design)) )) stop("All values must be finite and non-missing")
	design <- as.matrix(design)
	if(length(dispersion)<ngenes) dispersion <- rep(dispersion,length.out=ngenes)

	if(is.null(coef.start)) {
		start.method <- match.arg(start.method, c("null","y"))
		if(start.method=="null") N <- exp(offset)
	} else {
		coef.start <- as.matrix(coef.start)
	}
	
	offset <- t(expandAsMatrix(offset,dim(y)))

#	Check weights
	if(is.null(weights)) weights <- 1
	weights <- t(expandAsMatrix(weights,dim(y)))

# 	Initializing if desired. Note that lm.fit can fit in a vectorised manner, 
# 	where each column of the input matrix is a separate set of observations.
	if(is.null(coef.start)) {
		if(start.method=="y") {
			delta <- min(max(y), 1/6)
			y1 <- pmax(y, delta)
			fit <- lm.wfit(design, t(log(y1)) - offset, weights)
			beta <- fit$coefficients
			mu <- exp(beta + offset)
		} else {
			N <- expandAsMatrix(N,dim(y))
			w <- t(weights) * N/(1+dispersion*N)
			beta.mean <- log(rowSums(y*w/N)/rowSums(w))
			beta <- qr.coef(qr(design), matrix(beta.mean,nrow=nlibs,ncol=ngenes,byrow=TRUE))
			mu <- exp(design %*% beta + offset)
		}
	} else {
		beta <- t(coef.start)
		mu <- exp(design %*% beta + offset)
	}

# 	Checking arguments and calling the C++ method. We use transposed matrices so that each can be accessed from column-major storage in C++.
	if (!is.double(design)) storage.mode(design) <- "double"
	if (!is.double(dispersion)) storage.mode(dispersion) <- "double"
	if (!is.double(offset)) storage.mode(offset) <- "double"
	if (!is.double(weights)) storage.mode(weights) <- "double"
	if (!is.double(beta)) storage.mode(beta) <- "double"
	if (!is.double(mu)) storage.mode(mu) <- "double"
	output <- .Call("R_levenberg", nlibs, ngenes, design, t(y), dispersion, offset, weights, beta, mu, tol, maxit, PACKAGE="edgeR")

#	Check for error condition
	if (is.character(output)) { stop(output) }

#	Naming the output and returning it.  
	names(output) <- c("coefficients", "fitted.values", "deviance", "iter", "failed")
	output$coefficients <- t(output$coefficients)
	output$fitted.values <- t(output$fitted.values)
	colnames(output$coefficients) <- colnames(design)
	rownames(output$coefficients) <- rownames(y)
	dimnames(output$fitted.values) <- dimnames(y)
	output
}
