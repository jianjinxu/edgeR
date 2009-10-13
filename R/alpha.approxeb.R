alpha.approxeb<-function(object,verbose=TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to estimate a value for alpha, the weight for the weighted conditional log-likelihood, using approximate empirical Bayes methods
{
	if (!is(object,"DGEList"))
		stop("Currently supports DGEList objects")
	group<-object$samples$group
	levs.group<-levels(group)
	qA<-quantileAdjust(object,alpha=10,null.hypothesis=TRUE,verbose=verbose)  # alpha large to make common estimator
	d<-1/(1+qA$r[1])  # common delta
	scores<-0
	for(i in 1:length(levs.group)) {
		if (sum( group==levs.group[i] ) > 1) {
			#cat(levs.group[i],sum(group==levs.group[i]),"\n")
			scores<-scores+condLogLikDerDelta(qA$pseudo[,group==levs.group[i]],d,der=1,doSum=TRUE)
		}
	}
	exp.inf<-approx.expected.info(object,d,qA$pseudo)
	sigma2.0.est<-optimize(tau2.0.objective,c(0,500),info.g=exp.inf,score.g=scores)$min
	alpha<-1/(sigma2.0.est*nrow(object$counts)*mean(exp.inf))
	new("EBList",list(sigma2.0.est=sigma2.0.est,alpha=alpha,scores=scores,infos=exp.inf,quantileAdjusted=qA, common.dispersion.delta=d, common.dispersion=d/(1-d)))
}