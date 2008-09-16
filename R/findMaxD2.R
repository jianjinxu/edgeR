
findMaxD2<-function(x,alpha=0.5,grid=TRUE,tol=1e-05,n.iter=5,grid.length=200) {
  # this calculates delta as combined version of overall and individual
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]
  y2<-x$data[,x$group==g[2]]
  if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
    grid.vals<-seq(0.001,0.999,length=grid.length)
    l0<-condLogLikDerDelta(y1,grid.vals,der=0,doSum=TRUE)+condLogLikDerDelta(y2,grid.vals,der=0,doSum=TRUE)
    m0<-outer(rep(1,nrow(x$data)),l0)
    l0a<-l0+alpha*m0
    return(grid.vals[apply(l0a,1,which.max)])
  } else {  # do Newton Rhapson
    xprev<-findMaxD2(x,grid=TRUE,alpha=0,grid.length=20)
    mx<-tol+1; iter<-0
    while( mx > tol & iter < n.iter ) {
      iter<-iter+1
      l1<-condLogLikDerDelta(y1,xprev,der=1,grid=FALSE,doSum=TRUE)+condLogLikDerDelta(y2,xprev,der=1,grid=FALSE,doSum=TRUE)
      l1<-l1+sum(alpha*l1)
      l2<-condLogLikDerDelta(y1,xprev,der=2,grid=FALSE,doSum=TRUE)+condLogLikDerDelta(y2,xprev,der=2,grid=FALSE,doSum=TRUE)
      l2<-l1+sum(alpha*l2)
      x<-xprev-l1/l2
      mx<-max(abs(x-xprev))
      xprev<-x
    }
  }
  x
}


