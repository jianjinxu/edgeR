dispCoxReidSplineTrend <- function(y, design, offset=NULL, degree = 5, subset=1000, method.optim="Nelder-Mead", trace=0)
#	Estimate spline trend dispersion
#	Gordon Smyth, Yunshun Chen, Davis McCarthy
#	16 Dec 2010.  Last modified 24 Mar 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	method.optim <- match.arg(method.optim, c("Nelder-Mead", "BFGS"))
	
	abundance <- mglmOneGroup(y,offset=offset)
	i <- systematicSubset(subset, abundance)

#	Spline basis
	require("splines")
	p1 <- (1:(degree-1))/degree
	knots1 <- quantile(abundance,p=p1)
	r <- range(abundance)
	knots2 <- r[1]+p1*(r[2]-r[1])
	knots <- 0.3*knots1+0.7*knots2
	X <- cbind(1, ns(abundance, degree, knots=knots))
	
	fun <- function(par,y,design,offset,abundance,X) {
		eta <- X %*% par
		dispersion <- exp(eta - abundance)
		tryCatch(-sum(adjustedProfileLik(as.vector(dispersion),y,design,offset)),error=function(e) 1e10)
	}

	par0 <- rep(0,degree+1)
	par0[1] <- median(abundance[i]) + log(0.1)
	out <- optim(par0,fun,y=y[i,],design=design,offset=offset[i,],abundance=abundance[i],X=X[i,],control=list(trace=trace),method=method.optim)
	out$dispersion <- as.vector(exp(X %*% out$par - abundance))
	out$abundance <- abundance
	out
}
