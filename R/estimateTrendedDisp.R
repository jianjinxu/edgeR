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
	abundance.full <- log(object$conc$conc.common)
	group <- object$samples$group <- as.factor(object$samples$group)

	bindisp <- binCMLDispersion(object)

	if( method=="bin.spline" ) {
		require("splines")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(bindisp$abundance,p=p1)
		r <- range(bindisp$abundance)
		knots2 <- r[1]+p1*(r[2]-r[1])
		knots <- 0.3*knots1+0.7*knots2
		ind <- rep(NA, df+1)
		ind[1] <- which.min(bindisp$abundance)
		ind[df+1] <- which.max(bindisp$abundance)
		for(i in 2:df)
			ind[i] <- which.min(abs(knots[i-1]-bindisp$abundance))
		fit <- lm(dispersion ~ ns(abundance, df=df, knots=knots), data=bindisp)
		f <- splinefun(bindisp$abundance[ind], fit$fitted.value[ind], method="natural")
		dispersion <- f(abundance.full)
	}
	if( method=="bin.loess" ) {
		fit <- loessFit(bindisp$dispersion, bindisp$abundance, span=span)
		f <- approxfun(bindisp$abundance, fit$fitted, method="linear", rule=2)
		dispersion <- f(abundance.full)
	}

	object$trended.dispersion <- dispersion
	object
}

