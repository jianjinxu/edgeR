gof <- function( glmfit, pcutoff=0.1, adjust="holm", plot=FALSE, main="qq-plot of genewise goodness of fit", ...)
## Use LRT on deviance from a DGEGLM object
## to identify dispersion outlier genes
## Davis McCarthy, Gordon Smyth
## 23 Mar 2011. Last modified 23 May 2012.
{
	stopifnot( is(glmfit, "DGEGLM") )
	gof.stats <- glmfit$deviance
	gof.pvals <- pchisq(gof.stats, df=glmfit$df.residual, lower.tail=FALSE, log.p=FALSE)
	outlier <- p.adjust(gof.pvals, method=adjust) < pcutoff

	if(plot) {
		n <- length(gof.stats)
		z <- zscoreGamma(gof.stats,shape=glmfit$df.residual/2,scale=2)
		col <- rep("black",n)
		col[outlier] <- "blue"
		pch <- rep(1,n)
		pch[outlier] <- 16
		qqnorm(z,col=col,pch=pch,main=main,...)
		abline(0,1)
	}
	invisible(list(gof.statistics=gof.stats, gof.pvalues=gof.pvals, outlier=outlier, df=glmfit$df.residual[1]))
}


.gof2 <- function( y, design, dispersion, offset=NULL, pcutoff=0.1, fit=NULL, method="LR" )
    ## Calculate Deviance or Pearson goodness of fit statistics
    ## for the dispersion parameter and find dispersion outliers
    ## Davis McCarthy
    ## 8 Feb 2011. Last modified 13 Nov 2012.
{
    y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
    npar <- ncol(design)
    if( is.null(fit) )
        fit <- glmFit(y, design, dispersion, offset=offset, prior.count=0)
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



