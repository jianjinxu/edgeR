gof <- function( glmfit, pcutoff=0.1 )
    ## Use LRT on deviance from a DGEGLM object
    ## to identify dispersion outlier genes
    ## Davis McCarthy
    ## 23 Mar 2011. Last modified 23 Mar 2011.
{
    stopifnot( is(glmfit, "DGEGLM") )
    gof.stats <- glmfit$deviance
    gof.pvals <- pchisq(gof.stats, df=glmfit$df.residual, lower.tail=FALSE, log.p=FALSE)
    outlier <- p.adjust(gof.pvals, method="holm") < pcutoff

    new("list", list(gof.statistics=gof.stats, gof.pvalues=gof.pvals, outlier=outlier))
}


.gof2 <- function( y, design, dispersion, offset=NULL, pcutoff=0.1, fit=NULL, method="LR" )
    ## Calculate Deviance or Pearson goodness of fit statistics
    ## for the dispersion parameter and find dispersion outliers
    ## Davis McCarthy
    ## 8 Feb 2011. Last modified 23 Mar 2011.
{
    y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
    npar <- ncol(design)
    if( is.null(fit) )
        fit <- glmFit(y, design, dispersion, offset=offset)
    method <- match.arg(method, c("LR","Pearson"))

    if( method=="LR") {
        gof.stats <- fit$deviance
    }
    else {
        means <- fitted(fit)
        V <- means*(1+dispersion*means)
        gof.stats <- rowSums( (y-means)^2/V )
    }

    right <- gof.stats > ( nlibs - npar )
    gof.pvals <- rep(NA, ntags)
    gof.pvals[right] <- pchisq(gof.stats[right], df=(nlibs - npar), lower.tail=FALSE, log.p=FALSE)
    gof.pvals[!right] <- pchisq(gof.stats[!right], df=(nlibs - npar), lower.tail=TRUE, log.p=FALSE) 
    outlier <- p.adjust(gof.pvals, method="holm") < pcutoff

    new("list", list(gof.statistics=gof.stats, gof.pvalues=gof.pvals, outlier=outlier, right=right, fit=fit))
}



