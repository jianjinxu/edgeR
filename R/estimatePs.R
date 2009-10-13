estimatePs<-function(object,r,tol=1e-10,maxit=30) 
# Written by Mark Robinson, edited by Davis McCarthy, January 2009
# A function to estimate values of p (proportion) using iterative Newton's Method
# Expanded to multiple groups
# Let y be a list object, with elements corresponding to groups. Each element is a matrix for a particular group
{
	nrows<-nrow(object$counts)
	lib.size<-object$samples$lib.size
	group<-object$samples$group
	levs.group <- levels(group)
	onev<-rep(1,nrows)
	y<-splitIntoGroups(object)
	this.p.group<-matrix(0,nrow=nrows,ncol=length(levs.group), dimnames=list(NULL,levs.group))
	for(i in 1:length(levs.group)) {
		this.p.group[,i]<-rowMeans(y[[i]]/outer(onev,lib.size[group==levs.group[i]]))
	}
	this.p.com<-rowMeans(object$counts/outer(onev,lib.size))
	rsums<-matrix(0,nrow=nrows,ncol=length(levs.group))
	for(i in 1:length(levs.group)) {
		rsums[,i]<-rowSums(y[[i]])
	}
	min.val<-8.783496e-16
	for(i in 1:maxit) { # Newton-Raphson steps
		d1p.com<-logLikDerP(this.p.com,object$counts,lib.size,r,der=1)
		d1p<-matrix(0,nrow=nrows,ncol=length(levs.group))
		for(i in 1:length(levs.group)) {
			d1p[,i]<-logLikDerP(this.p.group[,i],y[[i]],lib.size[group==levs.group[i]],r,der=1)
			d1p[rsums[,i]==0,i]<-min.val
		}
		mx<-max(abs(c(d1p.com,d1p)))
		if( mx < tol ) { break }
		d2p<-matrix(0,nrow=nrows,ncol=length(levs.group))
		d2p.com<-logLikDerP(this.p.com,object$counts,lib.size,r,der=2)
		d2p<-matrix(0,nrow=nrows,ncol=length(levs.group))
		for(i in 1:length(levs.group)) {
			d2p[,i]<-logLikDerP(this.p.group[,i],y[[i]],lib.size[group==levs.group[i]],r,der=2)
			d2p[rsums[,i]==0,i]<-min.val
		}
		this.p.com<-this.p.com-d1p.com/d2p.com
		this.p.com[this.p.com<=0 | this.p.com>=1]<-1/(max(lib.size)*10)
		for(i in 1:length(levs.group)) {
			this.p.group[,i]<-this.p.group[,i]-d1p[,i]/d2p[,i]
			this.p.group[rsums[,i]==0,i]<-min.val
			this.p.group[(this.p.group[,i]<=0 | this.p.group[,i]>=1),i]<-1/(max(lib.size[group==levs.group[i]])*10)
		}
	}
	rownames(this.p.group)<-rownames(object$counts)
	return(list(conc.common=this.p.com,conc.group=this.p.group))
}