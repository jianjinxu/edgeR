dispCoxReidPowerTrend <- function(y, design, lib.size, abundance=NULL, trace=0)
#	Estimate trend dispersion=a*mean^b, log(dispersion)
#	Gordon Smyth
#	16 Dec 2010
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	offset <- matrix(log(lib.size),ntags,nlibs,byrow=TRUE)
	if(is.null(abundance)) abundance <- mglmOneGroup(y,offset=offset)
	
	fun <- function(par,y,design,offset,abundance) {
		dispersion <- exp(par[1]+par[2]*abundance)
		tryCatch(-sum(adjustedProfileLik(dispersion,y,design,offset)),error=function(e) 1e10)
	}

	par0 <- c(log(0.1),0)
	out <- optim(par0,fun,y=y,design=design,offset=offset,abundance=abundance,control=list(trace=trace))
	out$dispersion <- exp(out$par[1]+out$par[2]*abundance)
	out
}
