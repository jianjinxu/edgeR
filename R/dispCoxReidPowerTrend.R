dispCoxReidPowerTrend <- function(y, design, offset=NULL, subset=10000, method.optim="Nelder-Mead", trace=0)
#	Estimate trend dispersion=a*mean^b
#	Gordon Smyth, Davis McCarthy, Yunshun Chen
#	16 Dec 2010.  Last modified 14 Nov 2011.
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
		dispersion <- exp(par[1]+par[2]*abundance) + exp(par[3])
		tryCatch(-sum(adjustedProfileLik(dispersion,y,design,offset)),error=function(e) 1e10)
	}

	par0 <- c(log(0.1),0,-5)
	out <- optim(par0,fun,y=y.nonzero[i,],design=design,offset=offset.nonzero[i,],abundance=abundance.nonzero[i],control=list(trace=trace),method=method.optim)
	out$dispersion <- exp(out$par[1]+out$par[2]*abundance.full) + exp(out$par[3])
	out$abundance <- abundance.full
	out
}