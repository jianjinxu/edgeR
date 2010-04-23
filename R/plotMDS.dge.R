plotMDS.dge <- function (x, top=500, col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), ...) 
{
    if (is.matrix(x)){
          x <- DGEList(counts = x, group = c(rep("1", ncol(x))))
          colnames(x$counts) <- 1:dim(x$counts)[2]
    }
    if (!is(x, "DGEList"))
          stop("Currently supports DGEList object and matrix")
    mx <- as.matrix(x$counts)
    labels <- colnames(x$counts)
    nprobes <- nrow(mx)
    nsamples <- ncol(mx)
    if (ndim < 2) 
        stop("Need at least two dim.plot")
    if (nsamples < ndim) 
        stop("Too few samples")
    dd <- matrix(0, nrow = nsamples, ncol = nsamples)
    twd <- estimateTagwiseDisp(estimateCommonDisp(x), prior.n = 10, grid.length = 500)
    o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:min(nprobes,top)]
    subdata <- x$counts[o,]
    for (i in 2:(nsamples)) 
    	for (j in 1:(i - 1))
            dd[i, j] = sqrt(estimateCommonDisp(DGEList(counts=subdata[,c(i,j)], group=c("1","1")))$common.dispersion)
    a1 <- cmdscale(as.dist(dd), k = ndim)
    plot(a1[, dim.plot[1]], a1[, dim.plot[2]], type = "n", xlab = paste("Dimension",dim.plot[1]), ylab = paste("Dimension", dim.plot[2]), ... )
    text(a1[, dim.plot[1]], a1[, dim.plot[2]], labels = labels, col=col, cex=cex)
    return(invisible(dd))
}
