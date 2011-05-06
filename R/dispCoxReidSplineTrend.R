dispCoxReidSplineTrend <- function(y, design, offset=NULL, degree = 5, subset=1000, method.optim="Nelder-Mead", trace=0)
#	Estimate spline trend dispersion
#	Gordon Smyth, Yunshun Chen, Davis McCarthy
#	16 Dec 2010.  Last modified 4 May 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	method.optim <- match.arg(method.optim, c("Nelder-Mead", "BFGS"))

	all.zero <- rowSums(y)==0
    if( any(all.zero) )
        warning("Some rows of count matrix are all zero. These rows are ignored in dispersion calculation.")               
    abundance.full <- rep(NA, ntags)
	abundance.nonzero <- mglmOneGroup(y[!all.zero,],offset=offset[!all.zero,])
    abundance.full[all.zero] <- min(abundance.nonzero) - 0.5
    abundance.full[!all.zero] <- abundance.nonzero
	i <- systematicSubset(subset, abundance.nonzero)
    offset.nonzero <- offset[!all.zero,]
    y.nonzero <- y[!all.zero,]
    
#	Spline basis
	require("splines")
	p1 <- (1:(degree-1))/degree
	knots1 <- quantile(abundance.nonzero,p=p1)
	r <- range(abundance.nonzero)
	knots2 <- r[1]+p1*(r[2]-r[1])
	knots <- 0.3*knots1+0.7*knots2
	X <- cbind(1, ns(abundance.nonzero, degree, knots=knots))
	
	fun <- function(par,y,design,offset,abundance,X) {
		eta <- X %*% par
		dispersion <- exp(eta - abundance)
		tryCatch(-sum(adjustedProfileLik(as.vector(dispersion),y,design,offset)),error=function(e) 1e10)
	}

	par0 <- rep(0,degree+1)
	par0[1] <- median(abundance.nonzero[i]) + log(0.1)
	out <- optim(par0,fun,y=y.nonzero[i,],design=design,offset=offset.nonzero[i,],abundance=abundance.nonzero[i],X=X[i,],control=list(trace=trace),method=method.optim)
    disp <- rep(NA, ntags)
    disp.nonzero <-  as.vector(exp(X %*% out$par - abundance.nonzero))
    disp[all.zero] <- disp.nonzero[which.min(abundance.nonzero)]
    disp[!all.zero] <- disp.nonzero
    out$dispersion <- disp
	out$abundance <- abundance.full
	out
}
