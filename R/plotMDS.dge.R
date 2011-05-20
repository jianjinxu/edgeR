plotMDS.dge <- function (x, top=500, dd=NULL, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), ...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Last modified 20 May 2011
{
#   Check input
    if(is.matrix(x)) x <- DGEList(counts=x)
    if(!is(x,"DGEList")) stop("x must be a DGEList or a matrix")
    nprobes <- nrow(x)
    nsamples <- ncol(x)
    if(is.null(labels)) labels <- 1:nsamples
    if(ndim < 2) stop("dim.plot must be at least two")
    if(nsamples < ndim) stop("Too few samples")

    x$samples$group <- factor(rep.int(1,nsamples))
    if(!is.null(dd)){
        	if(dim(dd)[1]!= nsamples | dim(dd)[2]!= nsamples)
        		stop("The supplied distance matrix must have the same dimension as the object.")
    } else {
    	twd <- estimateTagwiseDisp(estimateCommonDisp(x), prior.n = 10, grid.length = 500)
    	o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:min(nprobes,top)]
    	subdata <- x$counts[o,]
    	dd <- matrix(0, nrow = nsamples, ncol = nsamples)
    	for (i in 2:(nsamples)) 
    	for (j in 1:(i - 1))
		dd[i, j] = sqrt(estimateCommonDisp(DGEList(counts=subdata[,c(i,j)]))$common.dispersion)
    }
    a1 <- cmdscale(as.dist(dd), k = ndim)
    plot(a1[, dim.plot[1]], a1[, dim.plot[2]], type = "n", xlab = paste("Dimension",dim.plot[1]), ylab = paste("Dimension", dim.plot[2]), ... )
    text(a1[, dim.plot[1]], a1[, dim.plot[2]], labels = labels, col=col, cex=cex)
    return(invisible(dd))
}
