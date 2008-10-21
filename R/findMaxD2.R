
findMaxD2<-function(x,alpha=0.5,grid=TRUE,tol=1e-05,n.iter=5,grid.length=200) {
  # this calculates delta as combined version of overall and individual
  g<-unique(x$group)
  y1<-x$data[,x$group==g[1]]
  y2<-x$data[,x$group==g[2]]
  onev<-rep(1,nrow(x$data))
  if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
    grid.vals<-seq(0.001,0.999,length=grid.length)
    l0<-condLogLikDerDelta(y1,grid.vals,der=0,doSum=FALSE)+condLogLikDerDelta(y2,grid.vals,der=0,doSum=FALSE)
    m0<-outer(onev,colSums(l0))
    l0a<-l0 + alpha*m0
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


.findMaxD2<-function(x,alpha=0.5,grid=TRUE,binary=!grid,tol=1e-05,n.iter=12,grid.length=200) {
  # this calculates delta as combined version of overall and individual
  g<-unique(x$group)
  y<-vector("list",length(g))
  onev<-rep(1,nrow(x$data))
  zerov<-rep(1,nrow(x$data))
  for(i in 1:length(y))
    y[[i]]<-x$data[,x$group==g[i]]
  if(grid & !binary) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
    grid.vals<-seq(0.001,0.999,length=grid.length)
    l0<-0
    for(i in 1:length(y))
	  if( sum(x$group==g[i]) > 1)
        l0<-l0+condLogLikDerDelta(y[[i]],grid.vals,grid=TRUE,der=0,doSum=FALSE)
    m0<-outer(onev,colSums(l0))
    l0a<-l0+alpha*m0
    return(grid.vals[apply(l0a,1,which.max)])
  } else if (binary & !grid) {  # do binary search
    grid<-outer(onev,c(0.001,0.999))
    for(j in 1:n.iter) {
	  left<-grid[,1]
	  right<-grid[,2]
      l0.left<-l0.right<-zerov
      for(i in 1:length(y)) {
	    if( sum(x$group==g[i]) > 1) {
          l0.left<-l0.left+condLogLikDerDelta(y[[i]],left,grid=FALSE,der=0,doSum=FALSE)
          l0.right<-l0.right+condLogLikDerDelta(y[[i]],right,grid=FALSE,der=0,doSum=FALSE)
		}
	  }
      l0a.left<-l0.left+alpha*sum(l0.left)
      l0a.right<-l0.right+alpha*sum(l0.right)
	  k<-l0a.left>l0a.right
	  mid<-rowSums(grid)/2
	  grid[k,]<-cbind(left[k],mid[k])
	  grid[!k,]<-cbind(mid[!k],right[!k])
	}
    return(grid.vals[apply(l0a,1,which.max)])
  } else if (!binary & !grid) {  # do Newton Rhapson
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
  } else { 
    stop("You have set grid=TRUE and binary=TRUE, pick one but not the other.")
  }
  x
}



