getOffset <- function(object)
#	Extract offset vector or matrix from data object and optional arguments
#	By default, offset is constructed from the lib.size and norm.factors
#	However offset supplied explicitly takes precedence
#	Arguments take precedence over corresponding values in object
#	Gordon Smyth
#	26 Jan 2011. Last modified 10 May 2011.
{
	offset <- object$offset
	lib.size <- object$samples$lib.size
	norm.factors <- object$samples$norm.factors
	
	if(!is.null(offset)) {
		return(offset)
	} else {		
		if(!is.null(norm.factors)) lib.size <- lib.size*norm.factors
		return(log(lib.size))
	}
}
