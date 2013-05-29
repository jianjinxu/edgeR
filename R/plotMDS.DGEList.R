plotMDS.DGEList <- function (x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), xlab=NULL, ylab=NULL, method="logFC", prior.count=2, gene.selection="pairwise", ...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 28 May 2013.
{
#	Check x
	x$counts <- as.matrix(x$counts)
	if(!all(is.finite(x$counts))) stop("Missing or infinite counts not allowed")

	nprobes <- nrow(x)
	nsamples <- ncol(x)
	if(nsamples < 3) stop("Need at least 3 columns of data")

#	Check value for labels
	if(is.null(labels)) labels <- 1:nsamples
	labels <- as.character(labels)

#	Check value for dim.plot
	if(nsamples < ndim) stop("Dimension to be plotted is greater than number of libraries")

#	Default method is to convert to moderated logCPM and call limma plotMDS
	method <- match.arg(method, c("logFC","bcv","BCV"))
	if(method=="logFC") {
		if(is.null(xlab)) xlab <- paste("Leading logFC dim",dim.plot[1])
		if(is.null(ylab)) ylab <- paste("Leading logFC dim",dim.plot[2])
		y <- cpm(x,log=TRUE,prior.count=prior.count)
		return(plotMDS(y,top=top,labels=labels,col=col,cex=cex,dim.plot=dim.plot,ndim=ndim,gene.selection=gene.selection,xlab=xlab,ylab=ylab,...))
	}

#	From here method="bcv"

	x$samples$group <- factor(rep.int(1,nsamples))

	cn <- colnames(x)
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))	

#	Check value for top
	if (top < nprobes) { 
		twd <- estimateTagwiseDisp(estimateCommonDisp(x), grid.length = 100) 
		o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:top]
		subdata <- x$counts[o,,drop=FALSE]
	} else {
		subdata<-x$counts
	}

	lib.size <- x$samples$lib.size * x$samples$norm.factors
	myFun <- function(delta, y, ...) sum(condLogLikDerDelta(y, delta, ...))

	for (i in 2:(nsamples)) {
		for (j in 1:(i - 1))  {
			mm <- subdata[,c(i,j)]
			rs5 <- rowSums(mm) > 5
			lib <- lib.size[c(i, j)]
			norm <- t(t(mm)/lib) * exp(mean(log(lib)))
			delta <- optimize(myFun, interval = c(0.0001,.99), tol = 0.000001, maximum = TRUE, y = norm[rs5,], der = 0)
			dd[i, j] = sqrt( delta$maximum / (1-delta$maximum) )
		}
	}

#	Securing against negative eigenvalues with non-Euclidian distance matrices.
	a1 <- cmdscale(as.dist(dd), k = ndim)
	ndiff <- ndim-ncol(a1)
	if (ndiff > 0) a1<-cbind(a1, matrix(runif(ndiff*nsamples, -1e-6, 1e-6), ncol=ndiff, nrow=nsamples))

	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	if(is.null(xlab)) xlab <- paste("BCV distance",dim.plot[1])
	if(is.null(ylab)) ylab <- paste("BCV distance",dim.plot[2])
	plotMDS(mds,labels=labels,col=col,cex=cex,xlab=xlab,ylab=ylab,...)
}
