estimateSmoothing<-function(object,verbose=TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to estimate a value for alpha, the weight for the weighted conditional log-likelihood, using approximate empirical Bayes methods
# Smoothing parameter is reported as "prior n", that is the weight given to the common conditional log-likelihood in terms of the number of observations that weight is equivalent to. E.g. having prior n = 10 means that the common conditional log-likelihood is given weight equivalent to 10 tag observations.
{
	if (!is(object,"DGEList"))
		stop("Currently supports DGEList objects")
	group<-object$samples$group
	levs.group<-levels(group)
	d<-object$common.dispersion/(object$common.dispersion + 1)  # common delta
	scores<-0
	for(i in 1:length(levs.group)) {
		if ( sum(group==levs.group[i]) > 1) {
			scores<-scores+condLogLikDerDelta(object$pseudo.alt[,group==levs.group[i]],d,der=1,grid=FALSE, doSum=TRUE)
		}
	}
	exp.inf<-approx.expected.info(object,d,object$pseudo.alt)
	names(exp.inf)<-rownames(object$counts)
	sigma2.0.est<-optimize(tau2.0.objective,c(0,500),info.g=exp.inf,score.g=scores)$min
	alpha<-1/(sigma2.0.est*nrow(object$counts)*mean(exp.inf))
	prior.n <- alpha*nrow(object$counts)
	prior.n
}