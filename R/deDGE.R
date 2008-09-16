
deDGE<-function(object,alpha=0.5,verbose=TRUE) {
  if (!is(object,"DGEList"))
    stop("Currently supports DGEList objects")
  g<-unique(object$group)
  if (length(g) > 2) { stop("Can only do 2 sample comparison here.") }
  k1<-object$group==g[1]; k2<-object$group==g[2]
  object$data<-as.matrix(object$data)
  if (verbose) cat("Calculating shrinkage overdispersion parameters.\n")
  qA<-quantileAdjust(object,alpha=alpha,verbose=verbose)
  r<-qA$r
  ps<-qA$ps
  pseudo<-qA$pseudo
  fisher<-exactTestNB(pseudo,object$group,qA$N*ps$p,r,verbose)
  rownames(pseudo)<-rownames(object$data)
  new("deDGEList",(list(ps=ps,r=r,pseudo=pseudo,M=qA$N,exact=fisher)))
}

