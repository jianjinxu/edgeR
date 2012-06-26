loessByCol <- function(y, x=NULL, weights=NULL, span=0.5, cell=0.01, method="loess")
#	Smooth columns of matrix by a simplified non-robust loess curve of degree 0.
#	Original Rcode version by Davis McCarthy, May 2010.
#	Original loess version by Yunshun Chen, 08 May 2012.
#	Modifications by Gordon Smyth.
#	Last modified 26 June 2012.
{
	y <- as.matrix(y)
	if(span<=0) return(y)
	if(span>1) span <- 1
	ntags <- nrow(y)
	if(is.null(x)) x <- 1:ntags

	method <- match.arg(method, c("Rcode","loess"))
	if(method=="loess") {
		if(is.null(weights)) weights <- rep(1,ntags)
		for(j in 1:ncol(y)) y[,j] <- unclass(stats:::simpleLoess(y=y[,j],x=x,weights=weights,degree=0,trace.hat="approximate",span=span,cell=cell))$fitted
		return(y)
	}

#	method="Rcode"
	if(!is.null(weights)) warning("weights ignored")
	radius <- floor(span*ntags/2)
	if(radius<1) return(y)
	nspan <- 2*radius+1
	if(ntags <= nspan) {
		left <- rep(1,ntags)
		right <- rep(ntags,ntags)
	} else {
		left <- pmax((1:ntags)-radius,1)
		right <- left+2*radius
		i <- right>ntags
		left[i] <- ntags-nspan+1
		right[i] <- ntags
	}
	o <- order(x)
	xo <- x[o]
	yo <- y[o,,drop=FALSE]
	for(i in 1:ntags) {
		window <- left[i]:right[i]
		deviation <- abs(xo[window]-xo[i])
		w <- deviation/max(deviation)
		w <- (1-w^3)^3
		w <- w/sum(w)
		y[i,] <- w %*% yo[window,]
	}
	oo <- order(o)
	y[oo,,drop=FALSE]
}
