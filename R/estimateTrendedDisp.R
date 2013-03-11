estimateTrendedDisp <- function(object, method="bin.spline", df=5, span=2/3)
# Yunshun Chen, Gordon Smyth.
# Created 2 Feb 2012, last modified on 4 Oct 2012.
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
		disp.bins[i] <- estimateCommonDisp(object[tagsinbin,])$common.dispersion
		logCPM.bins[i] <- mean(logCPM[tagsinbin])
	}

	if( method=="bin.spline" ) {
		require("splines")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(logCPM.bins,p=p1)
		r <- range(logCPM.bins)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		ind <- rep(NA, df+1)
		ind[1] <- which.min(logCPM.bins)
		ind[df+1] <- which.max(logCPM.bins)
		for(i in 2:df)
			ind[i] <- which.min(abs(knots[i-1]-logCPM.bins))
		fit <- lm(disp.bins ~ ns(logCPM.bins, df=df, knots=knots))
		f <- splinefun(logCPM.bins[ind], fit$fitted.value[ind], method="natural")
		dispersion <- f(logCPM)
	}
	if( method=="bin.loess" ) {
		fit <- loessFit(disp.bins, logCPM.bins, span=span)
		f <- approxfun(logCPM.bins, fit$fitted, method="linear", rule=2)
		dispersion <- f(logCPM)
	}

	object$trended.dispersion <- dispersion
	object
}

