gini <- function(x)
#	Gini diversity index for columns of a numeric matrix
#	Gordon Smyth
#	5 Feb 2016
{
	d <- dim(x)
	if(is.null(d)) d <- dim(x) <- c(length(x),1)
	for (j in 1:d[2]) x[,j] <- sort.int(x[,j],na.last=TRUE)
	(2*colSums((1:d[1])*x)/colSums(x)-d[1]-1)/d[1]
}
