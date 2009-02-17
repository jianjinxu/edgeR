# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to carry out the DE analysis of DGE data
deDGE<-function(object,alpha=500,doPoisson=FALSE,verbose=TRUE) {
	if (!is(object,"DGEList"))
		stop("Currently supports DGEList objects")
	g<-unique(object$group)
	#k1<-object$group==g[1]; k2<-object$group==g[2]
	object$data<-as.matrix(object$data)
	if(doPoisson) {
		if (verbose) cat("Quantile adjusting as Poisson.\n")
		qA<-quantileAdjust(object,r.init=1000,n.iter=1)
	} else {
		if (verbose) cat("Calculating shrinkage overdispersion parameters.\n")
		qA<-quantileAdjust(object,alpha=alpha,verbose=verbose)
	}
	rownames(qA$pseudo)<-rownames(object$data)
	new("deDGEList",(list(ps=qA$ps,r=qA$r,pseudo=qA$pseudo,group=object$group,M=qA$N)))
}

