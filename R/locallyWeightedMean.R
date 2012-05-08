locallyWeightedMean <- function(y, x, span=0.5, cell=0.01) {
## Function for faster tricube smoothing
## Written by Yunshun Chen, 08 May 2012.
	sl <- matrix(0, length(x), ncol(y))
	for(i in 1:ncol(y))
		sl[,i] <- unclass(stats:::simpleLoess(y[,i], x, weights=rep(1, length(x)), 
			degree=0, trace.hat="approximate", span=span, cell=cell))$fitted
	sl
}