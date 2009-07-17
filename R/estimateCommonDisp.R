# A function to estimate the common dispersion (using conditional maximum likelihood) for fixed counts (y), assuming library sizes are equal
# Calculated on the delta = phi/(1+phi) scale, returns dispersion on the phi and the delta scale
estimateCommonDisp <- function(object, grid=TRUE, tol=1e-05, n.iter=10, grid.length=200)
{
	nrows<-nrow(object$data)
	levs.group<-levels(object$group)
	y<-splitIntoGroups(object)
	if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		grid.vals<-seq(0.001,0.999,length=grid.length)
		l0<-0
		for(i in 1:length(y)) {
			l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0<-matrix(colSums(l0),nrow=1,ncol=grid.length)
		delta <- grid.vals[apply(m0,1,which.max)]
		return(list(dispersion=delta/(1-delta), dispersion.delta=delta)) # Returns estimate of common dispersion on phi scale and delta scale
	} else {  # do Newton Rhapson
		delta.prev<-estimateCommonDisp(object,grid=TRUE,grid.length=20) # Use rough grid search to get a starting estimate for the dispersion (delta scale)
		delta.prev<-delta.prev$dispersion.delta # Make the previous estimate of delta a numeric value on the delta scale
		if(!is.numeric(delta.prev)) {
			stop("delta.prev is not numeric")
		}
		mx<-tol+1; iter<-0
		while( mx > tol & iter < n.iter ) {
			iter<-iter+1
			l2<-l1<-0
			for(i in 1:length(y)) { # Find derivatives of conditional log-likelihood (genewise) by summing results for each group 
				l1<-l1+condLogLikDerDelta(y[[i]],delta.prev,der=1,grid=FALSE,doSum=TRUE)
				l2<-l2+condLogLikDerDelta(y[[i]],delta.prev,der=2,grid=FALSE,doSum=TRUE)
			}
			l1<-sum(l1) # Derivative of conditional log-likelihood for common dispersion obtained by summing over tags
			l2<-sum(l2) # Likewise, second derivative of conditional log-likelihood obtained by summing over tags
			delta<-delta.prev-l1/l2 # NR method to get new estimate of delta
			mx<-max(abs(delta-delta.prev))
			delta.prev<-delta
		}
		cat("Number of iterations was:  ",iter,"\n")
		if((delta < 0 | delta > 1)) {
			cat("Newton-Raphson estimate out of range, so fix using grid search\n")
			delta.fix<-estimateCommonDisp(object,grid=TRUE,grid.length=grid.length)
			delta <- delta.fix$dispersion.delta # Use dispersion estimate on delta scale from grid search as NR fails
		}
	}
	list(dispersion=delta/(1-delta), dispersion.delta=delta)
}