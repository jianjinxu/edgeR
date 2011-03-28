dispCoxReidPowerTrend <- function(y, design, offset=NULL, subset=1000, method.optim="Nelder-Mead", trace=0)
#	Estimate trend dispersion=a*mean^b
#	Gordon Smyth, Davis McCarthy
#	16 Dec 2010.  Last modified 08 Feb 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	method.optim <- match.arg(method.optim, c("Nelder-Mead", "BFGS"))
	
	abundance <- mglmOneGroup(y,offset=offset)
	i <- systematicSubset(subset, abundance)
	
	fun <- function(par,y,design,offset,abundance) {
		dispersion <- exp(par[1]+par[2]*abundance)
		tryCatch(-sum(adjustedProfileLik(dispersion,y,design,offset)),error=function(e) 1e10)
	}

	par0 <- c(log(0.1),0)
	out <- optim(par0,fun,y=y[i,],design=design,offset=offset[i,],abundance=abundance[i],control=list(trace=trace),method=method.optim)
	out$dispersion <- exp(out$par[1]+out$par[2]*abundance)
	out$abundance <- abundance
	out
}
