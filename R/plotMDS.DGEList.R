plotMDS.DGEList <- function (x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), xlab=paste("Dimension",dim.plot[1]), ylab=paste("Dimension",dim.plot[2]), ...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen and Gordon Smyth
#	23 May 2011.  Last modified 6 June 2011.
{
	library(edgeR)
#	Check input
	if(is.matrix(x)) x <- DGEList(counts=x)
    	if(!is(x,"DGEList")) stop("x must be a DGEList or a matrix")

#	Remove rows with missing or Inf values
	ok <- is.finite(x$counts)
	if(!all(ok)) x <- x[apply(ok,1,all),]
	if(is.null(labels)) labels<-1:dim(x)[2]

    	nprobes <- nrow(x)
    	nsamples <- ncol(x)

	if(is.null(labels)) labels <- 1:nsamples
	if(ndim < 2) stop("dim.plot must be at least two")
	if(nsamples < ndim) stop("Too few samples")

	x$samples$group <- factor(rep.int(1,nsamples))

	cn <- colnames(x)
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))	

	twd <- estimateTagwiseDisp(estimateCommonDisp(x), grid.length = 500)
    	o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:min(nprobes,top)]
    	subdata <- x$counts[o,]

	for (i in 2:(nsamples)) 
	for (j in 1:(i - 1))
		dd[i, j] = sqrt(estimateCommonDisp(DGEList(counts=subdata[,c(i,j)]))$common.dispersion)

	a1 <- cmdscale(as.dist(dd), k = ndim)
	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	plotMDS(mds,labels=labels,col=col,cex=cex,xlab=xlab,ylab=ylab,...)
}