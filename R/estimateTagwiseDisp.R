estimateTagwiseDisp <- function(object, prior.n=getPriorN(object), trend="movingave", span=0.5, prop.used=0.3, method="grid", grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE)
#  Tagwise dispersion using weighted conditional likelihood empirical Bayes.

#  Davis McCarthy, Mark Robinson, Yunshun Chen, Gordon Smyth.
#  Created 2009. Last modified 8 May 2012.
{
	if( !is(object,"DGEList") ) stop("object must be a DGEList")
	if( is.null(object$pseudo.alt) ) {
		message("Running estimateCommonDisp() on DGEList object before proceeding with estimateTagwiseDisp().")
		object <- estimateCommonDisp(object)
	}
	trend <- match.arg(trend,c("none","loess","movingave","tricube"))
	method <- match.arg(method,c("grid","optimize"))
	ntags <- nrow(object$counts)
	group <- object$samples$group <- as.factor(object$samples$group)
	y <- splitIntoGroups(list(counts=object$pseudo.alt,samples=object$samples))
	delta <- rep(0,ntags)

	if(method=="grid"){  # do spline interpolation
		if(verbose) message("Using interpolation to estimate tagwise dispersion. ")
		spline.pts <- seq(from=grid.range[1],to=grid.range[2],length=grid.length)
		spline.disp <- object$common.dispersion * 2^spline.pts
		grid.vals <- spline.disp/(1+spline.disp)
	
		l0 <- 0
		for(i in 1:length(y)) {
			l0 <- condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0 <- switch(trend,
			"loess" = ntags*locallyWeightedMean(l0, object$conc$conc.common, span=span),
			"movingave" = ntags*weightedComLikMA(object,l0,prop.used=span),
			"tricube" = ntags*weightedComLik(object,l0, prop.used=span),
			"none" = matrix(colSums(l0),ntags,grid.length,byrow=TRUE)
		)
		l0a <- l0 + prior.n/ntags*m0
		d <- rep(0,ntags)
		for(j in 1:ntags) d[j] <- maximizeInterpolant(spline.pts, l0a[j,])
		tagwise.dispersion <- object$common.dispersion * 2^d
	} else {	
		if(trend != "none") stop("optimize method doesn't allow for abundance-dispersion trend")
		if(verbose) message("Tagwise dispersion optimization begun, may be slow, progress reported every 100 tags")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0, doSum=FALSE)
			delta[tag] <- delta.this$maximum
			if(verbose) if(tag%%100==0) message("tag ",tag)
		}
		tagwise.dispersion <- delta/(1-delta)
	}
	if(verbose) cat("\n")
	
#	Output
	object$prior.n <- prior.n
	object$tagwise.dispersion <- tagwise.dispersion
	object
}

