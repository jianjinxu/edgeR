
deDGE<-function(object,alpha=500,doPoisson=FALSE,verbose=TRUE) {
  if (!is(object,"DGEList"))
    stop("Currently supports DGEList objects")
  g<-unique(object$group)
  if (length(g) > 2) { stop("Can only do 2 sample comparison here.") }
  k1<-object$group==g[1]; k2<-object$group==g[2]
  object$data<-as.matrix(object$data)
  if(doPoisson) {
    if (verbose) cat("Quantile adjusting as Poisson.\n")
    qA<-quantileAdjust(object,r.init=1000,n.iter=1)
  } else {
    if (verbose) cat("Calculating shrinkage overdispersion parameters.\n")
    qA<-quantileAdjust(object,alpha=alpha,verbose=verbose)
  }
  fisher<-exactTestNB(qA$pseudo,object$group,qA$N*qA$ps$p,qA$r,verbose)
  rownames(qA$pseudo)<-rownames(object$data)
  new("deDGEList",(list(ps=qA$ps,r=qA$r,pseudo=qA$pseudo,M=qA$N,exact=fisher)))
}

