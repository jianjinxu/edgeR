estimateTagwiseDisp <- function(object, prior.n=10, weighted.common=FALSE, prop.used=NULL, tol=1e-06, grid=TRUE, grid.length=200, verbose=TRUE)
    ## Written by Davis McCarthy, 2009. Last modified 11 June 2010.
    ## A function to estimate the common dispersion (using conditional maximum likelihood) for fixed counts (y), assuming library sizes are equal
    ## Must take equalized counts (pseudocounts), not raw counts
    ## Calculated on the delta = phi/(1+phi) scale, returns dispersion on the phi and the delta scale
    ## Now uses optimize instead of a grid search to estimate delta when not using NR methd - improves speed of function
{
	ntags<-nrow(object$counts)
	group<-object$samples$group<-as.factor(object$samples$group)
	levs.group<-levels(object$samples$group)
	y<-splitIntoGroups(new("DGEList",list(counts=object$pseudo.alt,samples=object$samples)))
	delta <- rep(0,ntags)
	onev<-rep(1,ntags)
	if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		if(verbose) cat("Using grid search to estimate tagwise dispersion. ")
		grid.vals<-seq(0.001,0.999,length.out=grid.length)
		l0<-0
		for(i in 1:length(y)) {
			l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
                if(weighted.common) {
                    if(is.null(prop.used))
                        prop.used <- 0.4
                    m0 <- ntags*weightedComLik(object,l0,prop.used=prop.used) # Weights sum to 1, so need to multiply by number of tags to give this the same weight overall as the regular common likelihood
                }
		else {
                    m0<-outer(onev,colSums(l0))
                }
                l0a<-l0 + prior.n/ntags*m0
		delta <- grid.vals[apply(l0a,1,which.max)]
	} else {	
		if(verbose) cat("Dispersion being estimated for tags (dot=1000 tags): ")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0, doSum=FALSE)
			delta[tag] <- delta.this$maximum
			if(verbose) {
				if(tag%%10==0) cat("'")
				if(tag%%1000==0) cat(".") 
			}
		}
	}
	if(verbose) cat("\n")
	tagwise.dispersion <- delta/(1-delta)
	new("DGEList",list(samples=object$samples, common.dispersion=object$common.dispersion, prior.n=prior.n, tagwise.dispersion=tagwise.dispersion, counts=object$counts, pseudo.alt=object$pseudo.alt, genes=object$genes, allZeros=object$allZeros, conc=object$conc, common.lib.size=object$common.lib.size))
}




.estimateTagwiseDisp <- function(object, prior.n=10, tol=1e-06, grid=TRUE, grid.length=1000, verbose=TRUE)
# A function to estimate the common dispersion (using conditional maximum likelihood) for fixed counts (y), assuming library sizes are equal
# Must take equalized counts (pseudocounts), not raw counts
# Calculated on the delta = phi/(1+phi) scale, returns dispersion on the phi and the delta scale
# Now uses optimize instead of a grid search to estimate delta when not using NR methd - improves speed of function
{
	ntags<-nrow(object$counts)
	levs.group<-levels(object$samples$group)
	y<-splitIntoGroups(object)
	delta <- rep(0,ntags)
	onev<-rep(1,ntags)
	if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		if(verbose) cat("Using grid search to estimate tagwise dispersion. ")
		grid.vals<-seq(0.001,0.999,length.out=grid.length)
		l0<-0
		for(i in 1:length(y)) {
			l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0<-outer(onev,colSums(l0))
		l0a<-l0 + prior.n/ntags*m0
		delta <- grid.vals[apply(l0a,1,which.max)]
	} else {	
		if(verbose) cat("Dispersion being estimated for tags (dot=1000 tags): ")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y, tag=tag, ntags=ntags, prior.n=prior.n, der=0, doSum=FALSE)
			delta[tag] <- delta.this$maximum
			if(verbose) {
				if(tag%%10==0) cat("'")
				if(tag%%1000==0) cat(".") 
			}
		}
	}
	if(verbose) cat("\n")
	list(dispersion=delta/(1-delta), dispersion.delta=delta) # Returns estimate of common dispersion on phi scale and delta scale
}
