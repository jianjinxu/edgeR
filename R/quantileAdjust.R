
quantileAdjust<-function(object,N=prod(object$lib.size)^(1/ncol(object$data)),alpha=0,null.hypothesis=FALSE,n.iter=5,r.init=NULL,tol=.001,verbose=TRUE) {
  # adjust data for estimation of phi (common dispersion only)
  if (is.null(r.init)) { 
    r.init<-1/findMaxD2(object,alpha=10)-1 
  }
  g<-unique(object$group)
  k1<-object$group==g[1]; k2<-object$group==g[2]
  y1<-matrix(object$data[,k1],ncol=sum(k1)); y2<-matrix(object$data[,k2],ncol=sum(k2))
  r<-r.init
  rprev<-r+1
  count<-0
  while( count < n.iter ) { 
    count<-count+1
    if (verbose) cat("[quantileAdjust] Iteration (dot=1000)",count,":")
    rprev<-r
    ps<-estimatePs(y1,y2,object$lib.size[k1],object$lib.size[k2],r)
    if (null.hypothesis==TRUE) {
      p1<-pnbinom(y1-1,size=r,mu=outer(ps$p,object$lib.size[k1]))+dnbinom(y1,size=r,mu=outer(ps$p,object$lib.size[k1]))/2
      p2<-pnbinom(y2-1,size=r,mu=outer(ps$p,object$lib.size[k2]))+dnbinom(y2,size=r,mu=outer(ps$p,object$lib.size[k2]))/2
      mu1<-outer(ps$p,rep(N,sum(k1)))
      mu2<-outer(ps$p,rep(N,sum(k2)))
    } else {
      p1<-pnbinom(y1-1,size=r,mu=outer(ps$p1,object$lib.size[k1]))+dnbinom(y1,size=r,mu=outer(ps$p1,object$lib.size[k1]))/2
      p2<-pnbinom(y2-1,size=r,mu=outer(ps$p2,object$lib.size[k2]))+dnbinom(y2,size=r,mu=outer(ps$p2,object$lib.size[k2]))/2
      mu1<-outer(ps$p1,rep(N,sum(k1)))
      mu2<-outer(ps$p2,rep(N,sum(k2)))
    }
    pseudo<-interpolateHelper(cbind(mu1,mu2),cbind(p1,p2),r,cbind(object$data[,k1],object$data[,k2]),verbose=verbose)
    pseudo[pseudo<0]<-0  # values less than zero for small samples seems to make this unstable
    d<-list(group=object$group,data=pseudo)
    r<-1/findMaxD2(d,alpha=alpha)-1
    if( max(abs(rprev-r)) < tol ) { break }
  }
  return(list(r=r,pseudo=pseudo,p=cbind(p1,p2),mu=cbind(mu1,mu2),ps=ps,N=N))
}

