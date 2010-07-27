estimateDispIter<-function(object,N=exp(mean(log(object$samples$lib.size))), prior.n=10,common.disp=FALSE, null.hypothesis=FALSE,n.iter=5,disp.init=NULL,tol=1e-6,verbose=TRUE) 
# Written by Davis McCarthy, September 2009, based on a function by Mark Robinson 
# A function to adjust counts for the estimation of common and tagwise dispersion(common dispersion only?)
# Returns estimate(s) of dispersion, pseudocounts, mean counts for tags, proportions for tags in samples and adjusted library size
{
	nrows<-nrow(object$counts)
	lib.size<-object$samples$lib.size
	group<-object$samples$group
	levs.group<-levels(group)
	y<-splitIntoGroups(object)
	if (is.null(disp.init)) { 
		q2q.pois <- equalizeLibSizes(object, disp=0, N=exp(mean(log(object$samples$lib.size))), null.hypothesis=null.hypothesis)
		pseudo.obj<-new("DGEList",list(counts=q2q.pois$pseudo, samples=object$samples))
		disp.init<- .estimateCommonDisp(pseudo.obj, tol=tol)
		disp.init<-rep(disp.init$dispersion, nrows)
	}
	if(length(disp.init)==1) {
		disp.init<-rep(disp.init, nrows)
	}
	disp<-disp.init
	disp.prev<-disp+1
	count<-0
	while( count < n.iter ) { 
		count<-count+1
		if (verbose) cat("Adjusting library sizes and estimating dispersion: iteration",count,"\n")
		disp.prev<-disp
		q2q.out<-equalizeLibSizes(object,disp=disp,null.hypothesis=null.hypothesis)  # Use Poisson distribution to equalize library sizes
		pseudo.obj<-new("DGEList",list(counts=q2q.out$pseudo, samples=object$samples))
		if(common.disp) {
			disp.obj <- .estimateCommonDisp(pseudo.obj, tol=tol)
			disp <- rep(disp.obj$dispersion, nrows)
		} else {
			disp.obj <- .estimateTagwiseDisp(pseudo.obj, prior.n=prior.n, tol=1e-06, grid=TRUE, grid.length=1000, verbose=TRUE)
			disp <- disp.obj$dispersion
		}
		if( max(abs(disp.prev-disp)) < tol ) { break }
	}
	return(list(dispersion=disp,pseudo=q2q.out$pseudo,mu=q2q.out$mu,conc=q2q.out$conc,N=N))
}
