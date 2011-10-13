estimateTagwiseDisp <- function(object, prior.n=getPriorN(object), trend="movingave", prop.used=0.3, method="grid", grid.length=200, tol=1e-06, verbose=TRUE)
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
	onev<-rep(1,ntags)
	if(method=="grid") {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		if(verbose) cat("Using grid search to estimate tagwise dispersion. ")
		grid.vals<-seq(0.001,0.999,length.out=grid.length)
		l0 <- 0
		for(i in 1:length(y)) {
			l0 <- condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0 <- switch(trend,
			# Weights sum to 1, so need to multiply by number of tags to give this the same weight overall as the regular common likelihood
			"movingave" = ntags*weightedComLikMA(object,l0,prop.used=prop.used),
			"tricube" = ntags*weightedComLik(object,l0,prop.used=prop.used),
			"none" = outer(onev,colSums(l0))
		)
		l0a<-l0 + prior.n/ntags*m0
		delta <- grid.vals[apply(l0a,1,which.max)]
	} else {	
		if(verbose) cat("Dispersion being estimated for tags (dot=1000 tags): ")
		if(trend != "none") warning("optimize method doesn't allow for abundance-dispersion trend")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0, doSum=FALSE)
			delta[tag] <- delta.this$maximum
			if(verbose) if(tag%%1000==0) cat(".")
		}
	}
	if(verbose) cat("\n")
	tagwise.dispersion <- delta/(1-delta)
	new("DGEList",list(samples=object$samples, common.dispersion=object$common.dispersion, prior.n=prior.n, tagwise.dispersion=tagwise.dispersion, counts=object$counts, pseudo.alt=object$pseudo.alt, genes=object$genes, all.zeros=object$all.zeros, conc=object$conc, common.lib.size=object$common.lib.size))
}

