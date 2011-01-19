dispCoxReidPowerTrend <- function(y, design, lib.size, offset=NULL, abundance=NULL, method.optim="Nelder-Mead", trace=0)
#	Estimate trend dispersion=a*mean^b
#	Gordon Smyth
#	16 Dec 2010.  Last modified 19 Jan 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	if(is.null(offset)) offset <- matrix(log(lib.size),ntags,nlibs,byrow=TRUE)
	method.optim <- match.args(method.optim, c("Nelder-Mead", "BFGS"))
	
	if(is.null(abundance)) abundance <- mglmOneGroup(y,offset=offset)
	
	fun <- function(par,y,design,offset,abundance) {
		dispersion <- exp(par[1]+par[2]*abundance)
		tryCatch(-sum(adjustedProfileLik(dispersion,y,design,offset)),error=function(e) 1e10)
	}

	par0 <- c(log(0.1),0)
	out <- optim(par0,fun,y=y,design=design,offset=offset,abundance=abundance,control=list(trace=trace),method=method.optim)
	out$dispersion <- exp(out$par[1]+out$par[2]*abundance)
	out$abundance <- abundance
	out
}
