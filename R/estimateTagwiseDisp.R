estimateTagwiseDisp <- function(object, prior.df=10, trend="movingave", span=NULL, method="grid", grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE)
#  Tagwise dispersion using weighted conditional likelihood empirical Bayes.
#  Davis McCarthy, Mark Robinson, Yunshun Chen, Gordon Smyth.
#  Created 2009. Last modified 11 March 2013.

#  Notes 3 July 2012:
#  - interpolating derivatives would be better than interpolating loglik values.
#  - share code with estimateGLMTagwiseDisp?
{
	if( !is(object,"DGEList") ) stop("object must be a DGEList")
	if(is.null(object$AveLogCPM)) object$AveLogCPM <- aveLogCPM(object)
	if( is.null(object$common.dispersion) ) {
		message("Running estimateCommonDisp() on DGEList object before proceeding with estimateTagwiseDisp().")
		object <- estimateCommonDisp(object)
	}
	trend <- match.arg(trend,c("none","loess","movingave","tricube"))
	if(trend=="tricube") trend="loess"
	method <- match.arg(method,c("grid","optimize"))
	ntags <- nrow(object$counts)
	group <- object$samples$group <- as.factor(object$samples$group)
	y <- splitIntoGroups(list(counts=object$pseudo.counts,samples=object$samples))
	delta <- rep(0,ntags)
	nlibs <- ncol(object$counts)
	ngroups <- length(unique(group))
	prior.n <- prior.df/(nlibs-ngroups)

	if(method=="grid"){  # do spline interpolation
		if(verbose) message("Using interpolation to estimate tagwise dispersion. ")
		spline.pts <- seq(from=grid.range[1],to=grid.range[2],length=grid.length)
		spline.disp <- object$common.dispersion * 2^spline.pts
		grid.vals <- spline.disp/(1+spline.disp)
	
		l0 <- matrix(0,ntags,grid.length)
		for(j in 1:grid.length) for(i in 1:length(y)) l0[,j] <- condLogLikDerDelta(y[[i]],grid.vals[j],der=0)+l0[,j]

		if(is.null(span)) if(trend=="movingave") span <- 0.3 else span <- 0.5
		m0 <- switch(trend,
 			"movingave" = {
 				o <- order(object$AveLogCPM)
 				oo <- order(o)
 				movingAverageByCol(l0[o,], width=floor(span*ntags))[oo,]
 			},
			"loess" = loessByCol(l0, object$AveLogCPM, span=span)$fitted.values,
			"none" = matrix(colMeans(l0),ntags,grid.length,byrow=TRUE)
		)
		l0a <- l0 + prior.n*m0
		d <- maximizeInterpolant(spline.pts, l0a)
		tagwise.dispersion <- object$common.dispersion * 2^d
	} else {	
		if(trend != "none") stop("optimize method doesn't allow for abundance-dispersion trend")
		if(verbose) message("Tagwise dispersion optimization begun, may be slow, progress reported every 100 tags")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0)
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

