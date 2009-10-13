# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to estimate the dispersion on the phi/(phi+1) scale
findMaxD2<-function(object,alpha=0.5,grid=TRUE,tol=1e-05,n.iter=10,grid.length=200) {
	# this calculates delta as combined version of overall and individual
	nrows<-nrow(object$counts)
	levs.group<-levels(object$samples$group)
	y<-splitIntoGroups(object)
	onev<-rep(1,nrows)
	if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		grid.vals<-seq(0.001,0.999,length=grid.length)
		l0<-0
		for(i in 1:length(y)) {
			l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0<-outer(onev,colSums(l0))
		l0a<-l0 + alpha*m0
		return(grid.vals[apply(l0a,1,which.max)])
	} else {  # do Newton Rhapson
		xprev<-findMaxD2(object,grid=TRUE,alpha=0,grid.length=20)
		if(!is.numeric(xprev)) {
			stop("xprev is not numeric")
		}
		if(length(xprev)!=nrows) {
			stop("delta is not of length nrow(y)")
		}
		mx<-tol+1; iter<-0
		while( mx > tol & iter < n.iter ) {
			iter<-iter+1
			l2<-l1<-0
			for(i in 1:length(y)) {
				l1<-l1+condLogLikDerDelta(y[[i]],xprev,der=1,grid=FALSE,doSum=TRUE)
				l2<-l2+condLogLikDerDelta(y[[i]],xprev,der=2,grid=FALSE,doSum=TRUE)
			}
			l1<-l1+sum(alpha*l1)
			l2<-l2+sum(alpha*l2)
			x<-xprev-l1/l2
			mx<-max(abs(x-xprev))
			xprev<-x
		}
		cat("Number of iterations was:  ",iter,"\n")
	}
	x
}

# A function to estimate the dispersion on the phi/(phi+1) scale
.findMaxD2b<-function(x,alpha=0.5,grid=FALSE,tol=1e-05,n.iter=10,grid.length=200) {
	# this calculates delta as combined version of overall and individual
	nrows<-nrow(x$data)
	y<-lapply(split(t(x$data),x$group), FUN=function(u) matrix(u,nrow=nrows,byrow=TRUE))
	onev<-rep(1,nrows)
	zerov<-rep(0,nrows)
	if(grid) {  # do a grid search, since some likelihoods may be monotone, not amenable to NR
		grid.vals<-seq(0.001,0.999,length=grid.length)
		l0<-0
		for(i in 1:length(y)) {
			l0<-condLogLikDerDelta(y[[i]],grid.vals,der=0,doSum=FALSE)+l0
		}
		m0<-outer(onev,colSums(l0))
		l0a<-l0 + alpha*m0
		return(grid.vals[apply(l0a,1,which.max)])
	} else {  # find maximums using optimize
		xprev<-mean(.findMaxD2b(x,grid=TRUE,alpha=0,grid.length=20))
		cat("One . for each call to weightLogLikDerDelta\n")
		x <- optimize(f=.weightLogLikDerDelta,y=y,der=0,alpha=alpha,doSum=TRUE,maximum=TRUE,interval=c(0,1))
		cat("\n")
	}
	x
}


.weightLogLikDerDelta<-function(x,y,der=1,alpha=0.5,doSum=TRUE) {
	cat(".")
	l<-0
	for(i in 1:length(y)) {
		l<-l+condLogLikDerDelta(y[[i]],delta=x,der=der,grid=FALSE,doSum=doSum)
	}
	wl<-l+sum(alpha*l)
	wl
}



# Use a binary search to ensure that there is a root in the interval
		#mx<-tol+1; iter<-0; rootfound<-FALSE
		#left<-zerov
		#right<-onev
		#mid<-rep(0.5,nrows)
		#while( !rootfound & iter < n.iter ) {
		#	iter<-iter+1
		#	l1right<-l1left<-l1mid<-0
		#	for(i in 1:length(y)) {
		#		l1left<-l1left+condLogLikDerDelta(y[[i]],left,der=1,grid=FALSE,doSum=TRUE)
		#		l1right<-l1right+condLogLikDerDelta(y[[i]],right,der=1,grid=FALSE,doSum=TRUE)
		#		l1mid<-l1mid+condLogLikDerDelta(y[[i]],mid,der=1,grid=FALSE,doSum=TRUE)
		#	}
		#	if( (l1left > 0) & (l1mid < 0) ) {
		#		rootfound<-TRUE
		#		xprev<-(left+mid)/2
		#	} else if( (l1mid > 0) & (l1right < 0) ) {
		#		rootfound<-TRUE
		#		xprev<-(right+mid)/2
		#	} else {
		#		
		#	}
		#}


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
    }
  } else { 
    stop("You have set grid=TRUE and binary=TRUE, pick one but not the other.")
  }
  x
}



