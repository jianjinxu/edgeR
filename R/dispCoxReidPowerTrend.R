dispCoxReidPowerTrend <- function(y, design, offset=NULL, subset=1000, method.optim="Nelder-Mead", trace=0)
#	Estimate trend dispersion=a*mean^b
#	Gordon Smyth, Davis McCarthy
#	16 Dec 2010.  Last modified 3 May 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	method.optim <- match.arg(method.optim, c("Nelder-Mead", "BFGS"))

    all.zero <- rowSums(y)==0
    if( any(all.zero) )
        warning("Some rows of count matrix are all zero. These rows are ignored in dispersion calculation.")               
    abundance.full <- rep(NA, ntags)
	abundance.nonzero <- mglmOneGroup(y[!all.zero,],offset=offset[!all.zero,])
    abundance.full[all.zero] <- min(abundance.nonzero) - 0.5
    abundance.full[!all.zero] <- abundance.nonzero
	i <- systematicSubset(subset, abundance.nonzero)
    offset.nonzero <- offset[!all.zero,]
    y.nonzero <- y[!all.zero,]
    
	fun <- function(par,y,design,offset,abundance) {
		dispersion <- exp(par[1]+par[2]*abundance)
		tryCatch(-sum(adjustedProfileLik(dispersion,y,design,offset)),error=function(e) 1e10)
	}

	par0 <- c(log(0.1),0)
	out <- optim(par0,fun,y=y.nonzero[i,],design=design,offset=offset.nonzero[i,],abundance=abundance.nonzero[i],control=list(trace=trace),method=method.optim)
	out$dispersion <- exp(out$par[1]+out$par[2]*abundance.full)
	out$abundance <- abundance.full
	out
}
