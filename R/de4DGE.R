de4DGE<-function(object,prior.n=10, disp.init=NULL, doPoisson=FALSE, useCommonDisp=TRUE, verbose=TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to carry out the DE analysis of DGE data
{
	if (!is(object,"DGEList"))
		stop("Currently supports DGEList objects")
	object$counts<-as.matrix(object$counts)
	if(doPoisson) {
		if (verbose) cat("Quantile adjusting as Poisson.\n")
		disp.out<-equalizeLibSizes(object,null.hypothesis=FALSE)
        dispersion <- "Not applicable for Poisson model"
	} else {
		if (verbose) cat("Calculating shrinkage overdispersion parameters.\n")
		disp.out<-estimateDispIter(object,common.disp=useCommonDisp,prior.n=prior.n,disp.init=disp.init,verbose=verbose)
		dispersion <- disp.out$dispersion
	}
	rownames(disp.out$pseudo)<-rownames(object$counts)
	colnames(disp.out$pseudo)<-paste("pseudo",colnames(object$counts),sep=".")
	new("de4DGEList",(list(conc=disp.out$conc,dispersion=dispersion,pseudo=disp.out$pseudo,group=object$samples$group,M=disp.out$N)))
}