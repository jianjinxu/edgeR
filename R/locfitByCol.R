locfitByCol <- function(y, x=NULL, weights=1, span=0.5, degree=0)
#	Gordon Smyth
#	20 Aug 2012.
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(is.null(x)) x <- 1:ntags
	if(span*ntags<2 || ntags<=1) return(y)
	suppressPackageStartupMessages(require(locfit))
	X <- cbind(1,x)
	for (j in 1:ncol(y)) y[,j] <- fitted(locfit.raw(X,y[,j],weights=weights, alpha=span,deg=degree))
	y
}
