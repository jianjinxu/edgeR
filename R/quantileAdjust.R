# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to adjust data for the estimation of phi (common dispersion only?)
quantileAdjust<-function(object,N=prod(object$lib.size)^(1/ncol(object$data)),alpha=0,null.hypothesis=FALSE,n.iter=5,r.init=NULL,tol=.001,verbose=TRUE) {
	# adjust data for estimation of phi (common dispersion only)
	if (is.null(r.init)) { 
		r.init<-1/findMaxD2(object,alpha=10)-1 
	}
	nrows<-nrow(object$data)
	lib.size<-object$lib.size
	group<-object$group
	levs.group<-levels(group)
	y<-splitIntoGroups(object)
	p<-matrix(0,nrow=nrows,ncol=ncol(object$data))
	mu<-matrix(0,nrow=nrows,ncol=ncol(object$data))
	r<-r.init
	rprev<-r+1
	count<-0
	while( count < n.iter ) { 
		count<-count+1
		if (verbose) cat("[quantileAdjust] Iteration (dot=1000)",count,":")
		rprev<-r
		ps<-estimatePs(object,r)
		if (null.hypothesis==TRUE) {
			for(i in 1:length(levs.group)) {
				p[,group==levs.group[i]]<-pnbinom(y[[i]]-1,size=r,mu=outer(ps$p.common,lib.size[group==levs.group[i]]))+dnbinom(y[[i]],size=r,mu=outer(ps$p.common,lib.size[group==levs.group[i]]))/2
				mu[,group==levs.group[i]]<-outer(ps$p.common,rep(N,sum(group==levs.group[i])))
			}
		} else {
			for(i in 1:length(levs.group)) {
				p[,group==levs.group[i]]<-pnbinom(y[[i]]-1,size=r,mu=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))+dnbinom(y[[i]],size=r,mu=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))/2
				mu[,group==levs.group[i]]<-outer(ps$p.group[,i],rep(N,sum(group==levs.group[i])))
			}
		}
		#pseudo<-interpolateHelper(cbind(mu1,mu2),cbind(p1,p2),r,cbind(object$data[,k1],object$data[,k2]),verbose=verbose)
		count.max<-apply(object$data,1,max)
		pseudo<-interpolateHelper(mu,p,r,count.max,verbose=verbose)
		pseudo[pseudo<0]<-0  # values less than zero for small samples seems to make this unstable
		d<-list(group=group,data=pseudo)
		r<-1/findMaxD2(d,alpha=alpha)-1
		if( max(abs(rprev-r)) < tol ) { break }
	}
	return(list(r=r,pseudo=pseudo,p=p,mu=mu,ps=ps,N=N))
}

