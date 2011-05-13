dispBinTrend <- function(y, design, offset=NULL, degree = 10, span=0.3, min.n=500, method.bin="CoxReid", method.trend="spline", trace=0, ...)
    ## Estimate dispersions with a trend by computing common dispersion Cox-Reid estimates of the dispersion in bins based on abundance and then fit a curve through the common dispersions using either natural cubic splines or loess. From the fitted curve, estimates of the dispersion for each tag/gene are obtained.
    ## Davis McCarthy, Gordon Smyth
    ## 10 Feb 2011.  Last modified 3 May 2011.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))

	method.bin <- match.arg(method.bin, c("CoxReid", "Pearson", "deviance"))
	method.trend <- match.arg(method.trend, c("spline", "loess"))

    abundance.full <- mglmOneGroup(y,offset=offset)
    abundance.full[is.infinite(abundance.full)] <- min(abundance.full[is.finite(abundance.full)], na.rm=TRUE) - 0.5

    bindisp <- binGLMDispersion( y, design, min.n=min.n, offset=offset, method=method.bin, ...) 

	if( method.trend=="spline" ) {
            ## Spline approach
            ## Spline basis
	    require("splines")
	    p1 <- (1:(degree-1))/degree
	    knots1 <- quantile(bindisp$abundance,p=p1)
	    r <- range(bindisp$abundance)
	    knots2 <- r[1]+p1*(r[2]-r[1])
	    knots <- 0.3*knots1+0.7*knots2
            fit <- lm(dispersion ~ ns(abundance, degree, knots=knots), data=bindisp)
	    dispersion <- predict(fit, data.frame(abundance=abundance.full))
	}
	if( method.trend=="loess" ) {
            ## Loess approach
            bindisp$abundance <- c(min(abundance.full), bindisp$abundance, max(abundance.full))
            bindisp$dispersion <- c(bindisp$dispersion[1], bindisp$dispersion, tail(bindisp$dispersion, n=1) )
            fit <- loess(dispersion ~ abundance, bindisp, span=span)
            dispersion <- predict(fit, data.frame(abundance=abundance.full), se=TRUE)$fit
	}
	neg <- dispersion < 0
    dispersion[neg] <- min(bindisp$dispersion)
    new("list", list(abundance=abundance.full, dispersion=dispersion, bin.abundance=bindisp$abundance, bin.dispersion=bindisp$dispersion))
}

