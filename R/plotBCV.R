plotBCV <- function( object , xlab="Abundance (log2 counts per million)", ylab="Biological coefficient of variation (BCV)", panel.first=grid(), ... ) 
    ## Make a plot of biological coefficient of variation against abundance (counts per million)
    ## Davis McCarthy
    ## Created 18 January 2012. Last modified 29 January 2012.
{
     if( !is(object, "DGEList") )
        stop("First argument must be a DGEList object.")
     if( is.null( object$tagwise.dispersion ) )
         stop("Tagwise dispersion estimates not found: object$tagwise.dispersion is NULL. Run either estimateTagwiseDisp() or estimateGLMTrendedDisp() then estimateGLMTagwiseDisp() before trying to plot BCV estimates.")
     if( is.null( object$abundance ) )
         object$abundance <- mglmOneGroup(as.matrix(object), dispersion=object$tagwise.dispersion)
     abund.cpm <- log2(exp(object$abundance)) + log2(1e6)
     plot(abund.cpm, sqrt(object$tagwise.dispersion), panel.first=panel.first, ylab=ylab, xlab=xlab )
     if( !is.null(object$trended.dispersion) ) {
         o <- order(abund.cpm)
         lines(abund.cpm[o], sqrt(object$trended.dispersion)[o], col="dodgerblue", lwd=3, lty=3)
     }
}


