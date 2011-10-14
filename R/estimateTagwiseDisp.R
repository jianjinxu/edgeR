estimateTagwiseDisp <- function(object, prior.n=getPriorN(object), trend="movingave", prop.used=0.3, method="grid", grid.length=200, tol=1e-06, verbose=FALSE)
# Tagwise dispersion using weighted conditional likelihood empirical Bayes.

# Davis McCarthy, Mark Robinson, Gordon Smyth.
# Created 2009. Last modified 14 Oct 2011.
{
	if( !is(object,"DGEList") ) stop("object must be a DGEList")
	if( is.null(object$pseudo.alt) ) {
		message("Running estimateCommonDisp() on DGEList object before proceeding with estimateTagwiseDisp().")
		object <- estimateCommonDisp(object)
	}
	trend <- match.arg(trend,c("none","movingave","tricube"))
	method <- match.arg(method,c("grid","optimize"))
	ntags <- nrow(object$counts)
	group <- object$samples$group <- as.factor(object$samples$group)
	y <- splitIntoGroups(list(counts=object$pseudo.alt,samples=object$samples))
	delta <- rep(0,ntags)
	if(method=="grid") {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		if(verbose) message("Using grid search to estimate tagwise dispersion. ")
		grid.vals<-seq(0.001,0.999,length.out=grid.length)
		l0 <- 0
		for(i in 1:length(y)) {
			l0 <- condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0 <- switch(trend,
			# Weights sum to 1, so need to multiply by number of tags to give this the same weight overall as the regular common likelihood
			"movingave" = ntags*weightedComLikMA(object,l0,prop.used=prop.used),
			"tricube" = ntags*weightedComLik(object,l0,prop.used=prop.used),
			"none" = matrix(colSums(l0),ntags,grid.length,byrow=TRUE)
		)
		l0a <- l0 + prior.n/ntags*m0
		delta <- grid.vals[apply(l0a,1,which.max)]
	} else {	
		if(trend != "none") stop("optimize method doesn't allow for abundance-dispersion trend")
		if(verbose) message("Tagwise dispersion optimization begun, may be slow, progress reported every 100 tags")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0, doSum=FALSE)
			delta[tag] <- delta.this$maximum
			if(verbose) if(tag%%100==0) message("tag ",tag)
		}
	}
	if(verbose) cat("\n")
	tagwise.dispersion <- delta/(1-delta)

#	Output
	object$prior.n <- prior.n
	object$tagwise.dispersion <- tagwise.dispersion
	object
}

