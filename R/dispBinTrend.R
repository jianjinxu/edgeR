dispBinTrend <- function(y, design, offset=NULL, df=5, span=2/3, min.n=500, method.bin="CoxReid", method.trend="spline", trace=0, ...)
#	Estimate common dispersion in bins based on abundance,
#	then fit a curve through the dispersions
#	Davis McCarthy, Gordon Smyth
#	Created 10 Feb 2011.  Last modified 20 Aug 2012.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))

	method.bin <- match.arg(method.bin, c("CoxReid", "Pearson", "deviance"))
	method.trend <- match.arg(method.trend, c("spline", "loess"))

	abundance.full <- mglmOneGroup(y,offset=offset,dispersion=0.01)
	abundance.full[is.infinite(abundance.full)] <- min(abundance.full[is.finite(abundance.full)], na.rm=TRUE) - 0.5

	bindisp <- binGLMDispersion( y, design, min.n=min.n, offset=offset, method=method.bin, ...) 
	nbins <- length(bindisp$dispersion)

	if(nbins < 7) {
		f <- approxfun(bindisp$abundance, sqrt(bindisp$dispersion), method="linear", rule=2)
		dispersion <- f(abundance.full)^2
		return(list(abundance=abundance.full, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion))
	}
	if( method.trend=="spline" ) {
		require("splines")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(bindisp$abundance,p=p1)
		r <- range(bindisp$abundance)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		ind <- rep(NA, df+1)
		ind[1] <- which.min(bindisp$abundance)
		ind[df+1] <- which.max(bindisp$abundance)
		for(i in 2:df) ind[i] <- which.min(abs(knots[i-1]-bindisp$abundance))
		fit <- lm.fit(ns(bindisp$abundance,df=df,knots=knots,intercept=TRUE), sqrt(bindisp$dispersion))
		f <- splinefun(bindisp$abundance[ind], fit$fitted.value[ind], method="natural")
		dispersion <- f(abundance.full)^2
	}
	if( method.trend=="loess" ) {
		fit <- loessFit(sqrt(bindisp$dispersion), bindisp$abundance, span=span)
		f <- approxfun(bindisp$abundance, fit$fitted, method="linear", rule=2)
		dispersion <- f(abundance.full)^2
	}
	minbindisp <- min(bindisp$dispersion)
	dispersion <- pmax(dispersion,minbindisp)
	list(abundance=abundance.full, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion)
}

