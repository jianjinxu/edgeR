estimateTrendedDisp <- function(object, method="bin.spline", df=5, span=2/3)
# Yunshun Chen, Gordon Smyth.
# Created 2 Feb 2012, last modified on 24 Sep 2015.
{
	if( !is(object,"DGEList") ) stop("object must be a DGEList")
	if( is.null(object$pseudo.counts) ) {
		message("Running estimateCommonDisp() on DGEList object before proceeding with estimateTrendedDisp().")
		object <- estimateCommonDisp(object)
	}

	ntags <- nrow(object$counts)	
	logCPM <- object$AveLogCPM
	if(is.null(logCPM)) logCPM <- aveLogCPM(object)
	
	nbins <- 50
	if(nbins>ntags) stop("nbins greater than number of rows of data")
	bins <- cutWithMinN(logCPM,intervals=nbins,min.n=floor(ntags/nbins))
	disp.bins <- logCPM.bins <- rep(NA,nbins)
	
	for(i in 1:nbins) {
		tagsinbin <- bins$group==i
		disp.bins[i] <- estimateCommonDisp(object[tagsinbin,],rowsum.filter=0)$common.dispersion
		logCPM.bins[i] <- mean(logCPM[tagsinbin])
	}

	if( method=="bin.spline" ) {
		if(!requireNamespace("splines",quietly=TRUE)) stop("splines required but is not available")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(logCPM.bins,p=p1)
		r <- range(logCPM.bins)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		basisbins <- splines::ns(logCPM.bins,df=df,knots=knots,intercept=TRUE)
		beta <- coefficients(lm.fit(basisbins, sqrt(disp.bins)))
		basisall <- predict(basisbins,newx=logCPM)
		dispersion <- drop(basisall %*% beta)^2
	}
	
	if( method=="bin.loess" ) {
		fit <- loessFit(sqrt(disp.bins), logCPM.bins, span=span, iterations=1)
		f <- approxfun(logCPM.bins, fit$fitted, rule=2)
		dispersion <- f(logCPM)^2
	}

	object$trended.dispersion <- dispersion
	object
}

