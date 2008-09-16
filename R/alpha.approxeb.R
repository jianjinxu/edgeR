

alpha.approxeb<-function(object,verbose=TRUE) {
  if (!is(object,"DGEList"))
    stop("Currently supports DGEList objects")
  g<-unique(object$group)
  k1<-object$group==g[1]; k2<-object$group==g[2]
  qA<-quantileAdjust(object,alpha=10,null.hypothesis=TRUE,verbose=verbose)  # alpha large to make common estimator
  d<-1/(1+qA$r[1])  # common delta
  scores<-condLogLikDerDelta(qA$pseudo[,k1],d,der=1,doSum=TRUE)+condLogLikDerDelta(qA$pseudo[,k2],d,der=1,doSum=TRUE)
  exp.inf<-approx.expected.info(object,d,qA)
  sigma2.0.est<-optimize(tau2.0.objective,c(0,500),info.g=exp.inf,score.g=scores)$min
  alpha<-1/(sigma2.0.est*nrow(object$data)*mean(exp.inf))
  new("EBList",list(sigma2.0.est=sigma2.0.est,alpha=alpha,scores=scores,infos=exp.inf,quantileAdjusted=qA))
}

