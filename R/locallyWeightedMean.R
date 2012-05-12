locallyWeightedMean <- function(y, x, weights=rep(1,nrow(y)), span=0.5, cell=0.01)
#	Smooth columns of matrix by non-robust loess curve of degree 0.
#	Written by Yunshun Chen, 08 May 2012.
{
	y <- as.matrix(y)
	for(j in 1:ncol(y))
		y[,j] <- unclass(stats:::simpleLoess(y=y[,j],x=x,weights=weights,degree=0,trace.hat="approximate",span=span,cell=cell))$fitted
	y
}
