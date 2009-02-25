# Written by Mark Robinson, edited by Davis McCarthy, January 2009
# A function to estimate values of p (proportion) using iterative Newton's Method
# Expanded to multiple groups
# Let y be a list object, with elements corresponding to groups. Each element is a matrix for a particular group
estimatePs<-function(object,r,tol=1e-10,maxit=30) {
	nrows<-nrow(object$data)
	lib.size<-object$lib.size
	group<-object$group
	#k<-unique(group)
	k <- levels(group)
	onev<-rep(1,nrows)
	y<-splitIntoGroups(object)
	this.p.group<-matrix(0,nrow=nrows,ncol=length(k), dimnames=list(NULL,k))
	for(i in 1:length(k)) {
		this.p.group[,i]<-rowMeans(y[[k[i]]]/outer(onev,lib.size[group==k[i]]))
	}
	this.p.com<-rowMeans(object$data/outer(onev,lib.size))
	rsums<-matrix(0,nrow=nrows,ncol=length(k))
	for(i in 1:length(k)) {
		rsums[,i]<-rowSums(y[[i]])
	}
	min.val<-8.783496e-16
	for(i in 1:maxit) { # do 10 Newton method steps
		d1p.com<-logLikDerP(this.p.com,object$data,lib.size,r,der=1)
		d1p<-matrix(0,nrow=nrows,ncol=length(k))
		for(i in 1:length(k)) {
			d1p[,i]<-logLikDerP(this.p.group[,i],y[[i]],lib.size[group==k[i]],r,der=1)
			d1p[rsums[,i]==0,i]<-min.val
		}
		mx<-max(abs(c(d1p.com,d1p)))
		if( mx < tol ) { break }
		d2p<-matrix(0,nrow=nrows,ncol=length(k))
		d2p.com<-logLikDerP(this.p.com,object$data,lib.size,r,der=2)
		d2p<-matrix(0,nrow=nrows,ncol=length(k))
		for(i in 1:length(k)) {
			d2p[,i]<-logLikDerP(this.p.group[,i],y[[i]],lib.size[group==k[i]],r,der=2)
			d2p[rsums[,i]==0,i]<-min.val
		}
		this.p.com<-this.p.com-d1p.com/d2p.com
		this.p.com[this.p.com<=0 | this.p.com>=1]<-1/(max(lib.size)*10)
		for(i in 1:length(k)) {
			this.p.group[,i]<-this.p.group[,i]-d1p[,i]/d2p[,i]
			this.p.group[rsums[,i]==0,i]<-min.val
			this.p.group[(this.p.group[,i]<=0 | this.p.group[,i]>=1),i]<-1/(max(lib.size[group==k[i]])*10)
		}
	}
	return(list(p.common=this.p.com,p.group=this.p.group))
}



# Written by Mark Robinson
# A function to estimate values of p (proportion) using iterative Newton's Method
.estimatePs<-function(y1,y2,lib.size1,lib.size2,r,tol=1e-10,maxit=30) {
	if (is.vector(y1) | is.vector(y2)) {
		y1<-matrix(y1,nrow=1)
		y2<-matrix(y2,nrow=1)
	}
	y<-cbind(y1,y2)
	onev<-rep(1,nrow(y1))
	lib.size<-c(lib.size1,lib.size2)
	this.p1<-rowMeans(y1/outer(onev,lib.size1))
	this.p2<-rowMeans(y2/outer(onev,lib.size2))
	this.p<-rowMeans(y/outer(onev,lib.size))
	a<-rowSums(y1)
	b<-rowSums(y2)
	min.val<-8.783496e-16
	for(i in 1:maxit) { # do 10 Newton method steps
		d1p<-logLikDerP(this.p,y,lib.size,r,der=1)
		d1p1<-logLikDerP(this.p1,y1,lib.size1,r,der=1); d1p1[a==0]<-min.val
		d1p2<-logLikDerP(this.p2,y2,lib.size2,r,der=1); d1p2[b==0]<-min.val
		mx<-max(abs(c(d1p,d1p1,d1p2)))
		if( mx < tol ) { break }
		d2p<-logLikDerP(this.p,y,lib.size,r,der=2)
		d2p1<-logLikDerP(this.p1,y1,lib.size1,r,der=2)
		d2p2<-logLikDerP(this.p2,y2,lib.size2,r,der=2)
		this.p<-this.p-d1p/d2p
		this.p1<-this.p1-d1p1/d2p1
		this.p1[a==0]<-min.val
		this.p2<-this.p2-d1p2/d2p2
		this.p2[b==0]<-min.val
		this.p[this.p<=0 | this.p>=1]<-1/(max(c(lib.size1,lib.size2))*10)
		this.p1[this.p1<=0 | this.p1>=1]<-1/(max(lib.size1)*10)
		this.p2[this.p2<=0 | this.p2>=1]<-1/(max(lib.size2)*10)
	}
	return(list(p=this.p,p1=this.p1,p2=this.p2))
}



