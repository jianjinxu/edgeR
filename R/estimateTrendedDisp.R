estimateTrendedDisp <- function(object, method="bin.spline", df=5, span=2/3)
# Yunshun Chen, Gordon Smyth.
# Created 2 Feb 2012.
{
	if( !is(object,"DGEList") ) stop("object must be a DGEList")
	if( is.null(object$pseudo.alt) ) {
		message("Running estimateCommonDisp() on DGEList object before proceeding with estimateTrendedDisp().")
		object <- estimateCommonDisp(object)
	}

	ntags <- nrow(object$counts)	
	logCPM <- object$logCPM

	bindisp <- binCMLDispersion(object)

	if( method=="bin.spline" ) {
		require("splines")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(bindisp$logCPM.bins,p=p1)
		r <- range(bindisp$logCPM.bins)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		ind <- rep(NA, df+1)
		ind[1] <- which.min(bindisp$logCPM.bins)
		ind[df+1] <- which.max(bindisp$logCPM.bins)
		for(i in 2:df)
			ind[i] <- which.min(abs(knots[i-1]-bindisp$logCPM.bins))
		fit <- lm(dispersion.bins ~ ns(logCPM.bins, df=df, knots=knots), data=bindisp)
		f <- splinefun(bindisp$logCPM.bins[ind], fit$fitted.value[ind], method="natural")
		dispersion <- f(logCPM)
	}
	if( method=="bin.loess" ) {
		fit <- loessFit(bindisp$dispersion.bins, bindisp$logCPM.bins, span=span)
		f <- approxfun(bindisp$logCPM.bins, fit$fitted, method="linear", rule=2)
		dispersion <- f(logCPM)
	}

	object$trended.dispersion <- dispersion
	object
}

