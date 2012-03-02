## DGLMSTDRESID.R

dglmStdResid <- function(y, design, dispersion=0, offset=0, nbins=100, make.plot=TRUE, xlab="Mean", ylab="Ave. binned standardized residual", ... ) {
    ## Function to bin DGE data based on fitted values for the abundance and compute and plot the average of the standardized residuals from a Poisson model fit in each bin against the average abundance for each bin. Allows us to investigate the mean-variance relationship in the data and compute a variance function for the negative binomial model.
    ## Davis McCarthy
    ## Created 9 November 2010. Last modified 1 March  2012.
	ngenes <- nrow(y)
	nlibs <- ncol(y)
    if( length(offset)!=nlibs & length(offset)!=1 & length(offset)!=length(y) )
		stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
	else
		offset <- matrix(offset, nrow=ngenes, ncol=nlibs, byrow=TRUE)
	fit <- mglmLS(y, design=design, dispersion=0, offset=offset)
    means <- as.vector(fit$fitted)
    std.resid <- nlibs * ( as.vector(y) - means )^2 / ( nlibs - ncol(design) ) # Obtain an approximate value for the standardized residual: denominator is (n - p) / n instead of the usual (1 - leverage)
    n <- length(means)
    means.quantiles <- quantile(means, probs=seq(0,1,length=nbins+1))
    means.quantiles[1] <- 0
    if(any(duplicated(means.quantiles)))
        stop("Duplicated quantiles for the means, so cannot produce bins. Try altering nbins.")
    else
        f <- cut(means,breaks=means.quantiles)
    bins <- split(1:n,f)
    std.resid.bins <- means.bins <- list()
    for(i in 1:nbins){
        means.bins[[i]] <- means[bins[[i]]]
        std.resid.bins[[i]] <- std.resid[bins[[i]]]
    }
    ave.means <- sapply(means.bins, mean)
    ave.std.resid <- sapply(std.resid.bins, mean)
    out <- list(ave.means=ave.means, ave.std.resid=ave.std.resid, bin.means=means.bins, bin.std.resid=std.resid.bins, means=means, standardized.residuals=std.resid, bins=bins, nbins=nbins, ngenes=ngenes, nlibs=nlibs)
    out$dispersion.estimate <- getDispersions(out)
    if(make.plot)
        plot(ave.means, ave.std.resid, pch="x", col="darkgreen", cex=1.5, log="xy", xlab=xlab, ylab=ylab, plot.first=grid(), ...)
    return( invisible( out ) )
}


getDispersions <- function(binned.object) {
    ## Estimate the dispersion parameter for each DGE observation from the variance function calculated by binning the standardized residuals from the Poisson GLM based on the estimated (fitted) mean for the observation. Operates on the output of binStdResidPois
    ## Davis McCarthy
    ## Created 9 November 2010. Last modified 9 November 2010.
    dispersion <- rep(NA, length=length(binned.object$means))
    bin.dispersion <- ( binned.object$ave.std.resid - binned.object$ave.means) / binned.object$ave.means^2
    bin.dispersion.used <- bin.dispersion
    whichbin <- 1:binned.object$nbins
    for( i in 1:binned.object$nbins) {
        if(bin.dispersion[i] > 0)
            dispersion[binned.object$bins[[i]]] <- bin.dispersion.used[i] <- bin.dispersion[i]
        else {
            next.ok <- min( whichbin[ bin.dispersion > 0 & whichbin > i ] )
            dispersion[binned.object$bins[[i]]] <- bin.dispersion.used[i] <- bin.dispersion[ next.ok ]
        }
    }
    dispersion <- matrix(dispersion, nrow=binned.object$ngenes, ncol=binned.object$nlibs)
    list(bin.dispersion=bin.dispersion, bin.dispersion.used=bin.dispersion.used, dispersion=dispersion)
}


    




