plotMDS.DGEList <- function (x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), xlab=paste("Dimension",dim.plot[1]), ylab=paste("Dimension",dim.plot[2]), ...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 23 July 2012.
{
#	Remove rows with missing or Inf values
	ok <- is.finite(x$counts)
	if(!all(ok)) x <- x[rowSums(ok)>0,]
	nprobes <- nrow(x)
	nsamples <- ncol(x)

#	Check value for top
	top <- min(top,nprobes)

#	Check value for labels
	if(is.null(labels)) labels <- 1:nsamples

#	Check value for dim.plot
	if(nsamples < ndim) stop("Dimension to be plotted is greater than number of libraries")

	x$samples$group <- factor(rep.int(1,nsamples))

	cn <- colnames(x)
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))	

	twd <- estimateTagwiseDisp(estimateCommonDisp(x), grid.length = 100) # <- unnecessary with spline interpolation?
	o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:min(nprobes,top)]
	subdata <- x$counts[o,]

	gm <- function(x) exp( mean(log(x)) )
	myFun <- function(delta, y, ...) sum(condLogLikDerDelta(y, delta, ...))

	for (i in 2:(nsamples)) {
		for (j in 1:(i - 1))  {
			mm <- subdata[,c(i,j)]
			rs5 <- rowSums(mm) > 5
			norm <- t( t(mm) / colSums(mm) )*gm(colSums(mm))
			delta <- optimize(myFun, interval = c(0.0001,.99), tol = 0.000001, maximum = TRUE, y = norm[rs5,], der = 0)
			dd[i, j] = sqrt( delta$maximum / (1-delta$maximum) )
		}
	}

#	Securing against negative eigenvalues with non-Euclidian distance matrices.
	a1 <- cmdscale(as.dist(dd), k = ndim)
	ndiff<-ndim-ncol(a1)
	if (ndiff > 0) a1<-cbind(a1, matrix(runif(ndiff*nsamples, -1e-6, 1e-6), ncol=ndiff, nrow=nsamples))

	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	plotMDS(mds,labels=labels,col=col,cex=cex,xlab=xlab,ylab=ylab,...)
}
