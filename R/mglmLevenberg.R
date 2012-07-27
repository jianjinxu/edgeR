mglmLevenberg <- function(y, design, dispersion=0, offset=0, coef.start=NULL, start.method="null", tol=1e-06, maxit=200)
#	Fit genewise negative binomial glms with log-link
#	using Levenberg damping to ensure secure convergence

#	R version by Gordon Smyth and Yunshun Chen
#	C++ version by Aaron Lun
#	Created 3 March 2011.  Last modified 11 July 2012
{
#	Check arguments
	y <- as.matrix(y)
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
		start <- as.matrix(start)
	}
	offset <- expandAsMatrix(offset,dim(y))

	# Initializing if desired. Note that lm.fit can fit in a vectorised manner, 
	# where each column of the input matrix is a separate set of observations.
	if(is.null(coef.start)) {
		if(start.method=="y") {
			delta <- min(max(y), 1/6)
			y1 <- pmax(y, delta)
			fit <- lm.fit(design, t(log(y1) - offset))
			beta <- t(fit$coefficients)
			mu <- exp(t(fit$fitted.values) + offset)
		} else {
			N <- expandAsMatrix(N,dim(y))
			beta.mean <- log(.rowMeans(y/N,ngenes,nlibs))
			beta <- qr.coef(qr(design), matrix(beta.mean,nrow=nlibs,ncol=ngenes,byrow=TRUE))
			mu <- exp(t(design %*% beta) + offset)
			beta <- t(beta)
		}
	} else {
		beta <- coef.start
		mu <- exp(beta %*% t(design) + offset)
	}

	# Calling the C++ method.
	output <- .Call("mglm_levenberg", nlibs, ngenes, design, y, dispersion, offset, beta, mu, tol, maxit, PACKAGE="edgeR")

	# Naming the output and returning it.  
	names(output) <- c("coefficients", "fitted.values", "deviance", "iter", "failed")
	colnames(output$coefficients) <- colnames(design)
	rownames(output$coefficients) <- rownames(y)
	dimnames(output$fitted.values) <- dimnames(y)
	output
}
