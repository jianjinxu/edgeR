plotBCV <- function(object, xlab="logCPM", ylab="Biological coefficient of variation", pch=16, cex=0.2, ...)
#	Plot of biological coefficient of variation against abundance (log CPM)
#	Davis McCarthy, Yunshun Chen, Gordon Smyth.
#	Created 18 January 2012.  Last modified 2 Nov 2012.
{
	if(!is(object,"DGEList")) stop("object must be a DGEList object.")
	if(is.null(object$tagwise.dispersion)) stop("Tagwise dispersion estimates not found. Run either estimateTagwiseDisp() or estimateGLMTrendedDisp() then estimateGLMTagwiseDisp() before trying to plot BCV estimates.")
#	if(is.null(object$abundance)) object$abundance <- mglmOneGroup(object$counts, offset=getOffset(object), dispersion=object$tagwise.dispersion)
	logCPM <- object$logCPM
	if(is.null(logCPM)) {
		abundance <- mglmOneGroup(object$counts,offset=getOffset(object),dispersion=0.05)
		logCPM <- log1p(exp(abundance+log(1e6)))/log(2)
	}
	plot(logCPM, sqrt(object$tagwise.dispersion), ylab=ylab, xlab=xlab, pch=pch, cex=cex, ...)
	if(!is.null(object$common.dispersion)) abline(h=sqrt(object$common.dispersion), col="red", lwd=2)	
	if(!is.null(object$trended.dispersion)) {
		o <- order(logCPM)
		lines(logCPM[o], sqrt(object$trended.dispersion)[o], col="blue", lwd=2)
		legend("topright",legend=c("Trend","Common"),lwd=2,col=c("blue","red"))
	} else {
		legend("topright",legend="Common",lwd=2,col="red")
	}
	invisible()
}
