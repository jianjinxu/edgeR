dispBinTrend <- function(y, design=NULL, offset=NULL, df=5, span=2/3, min.n=400, method.bin="CoxReid", method.trend="spline", trace=0, abundance=NULL, ...)
#	Estimate common dispersion in bins based on abundance,
#	then fit a curve through the dispersions
#	Davis McCarthy, Gordon Smyth
#	Created 10 Feb 2011.  Last modified 20 Aug 2012.
{
#	Check y
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,nlibs,1)
	} else {
		design <- as.matrix(design)
	}

#	Check offset
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))

	method.bin <- match.arg(method.bin, c("CoxReid", "Pearson", "deviance"))
	method.trend <- match.arg(method.trend, c("spline", "loess"))

	if(is.null(abundance)) {
		abundance <- mglmOneGroup(y,offset=offset,dispersion=0.02)
		abundance <- log2(exp(abundance+log(1e6))+0.5)
	}
	
	bindisp <- binGLMDispersion(y, design, min.n=min.n, offset=offset, method=method.bin, abundance=abundance, ...) 

	nbins <- length(bindisp$dispersion)
	if(nbins==1) {
		dispersion <- rep(bindisp$dispersion,ntags)
		return(list(abundance=abundance, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion))
	}
	logd <- log(bindisp$dispersion+0.01)
#	if(nbins==2) {
#		logd <- log(bindisp$dispersion+0.01)
#		a <- sum(logd)/2
#		b <- diff(logd)/diff(bindisp$abundance)
#		x <- abundance-mean(bindisp$abundance)
#		dispersion <- a+b*x
#		dispersion <- exp(dispersion)-0.01
#		return(list(abundance=abundance, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion))
#	}

	if(nbins<7) {
		f <- approxfun(bindisp$abundance,sqrt(bindisp$dispersion),rule=2)
		dispersion <- f(abundance)^2
#		dispersion <- sqrt(dispersion)-0.01
		return(list(abundance=abundance, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion))
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
		dispersion <- f(abundance)^2
	}

	if( method.trend=="loess" ) {
		fit <- loessFit(sqrt(bindisp$dispersion), bindisp$abundance, span=span)
		f <- approxfun(bindisp$abundance, fit$fitted, rule=2)
		dispersion <- f(abundance)^2
	}

	minbindisp <- min(bindisp$dispersion)
#	dispersion <- pmax(dispersion,minbindisp)
	list(abundance=abundance, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion)
}

