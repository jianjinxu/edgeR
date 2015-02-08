dropEmptyLevels <- function(x)
#	Drop levels of a factor that don't occur
#	Gordon Smyth
#	Created 25 March 2012.  Last modified 8 Feb 2015.
{
	if(is.factor(x)) {
		v <- as.integer(x)
		i <- which(tabulate(v)>0L)
		if(length(i) < nlevels(x)) {
			levels(v) <- levels(x)[i]
			class(v) <- class(x)
			return(v)
		} else {
			return(x)
		}
	} else {
		return(factor(x))
	}
}
