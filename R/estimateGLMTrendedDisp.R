#  Last modified 23 April 2012

estimateGLMTrendedDisp <- function(y, design=NULL, offset=NULL, method="auto", ...) 
UseMethod("estimateGLMTrendedDisp")

estimateGLMTrendedDisp.DGEList <- function(y, design=NULL, offset=NULL, method="auto", ...)
{
    if( is.null(offset) )
        offset <- getOffset(y)
	d <- estimateGLMTrendedDisp(y=y$counts, design=design, offset=offset, method=method, ...)
	y$trended.dispersion <- d$dispersion
	y$abundance <- d$abundance
	if( !is.null(d$bin.dispersion) ) y$bin.dispersion <- d$bin.dispersion
	if( !is.null(d$bin.abundance) ) y$bin.abundance <- d$bin.abundance
	y
}

estimateGLMTrendedDisp.default <- function(y, design=NULL, offset=NULL, method="auto", ...)
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: cannot estimate dispersion")
		return(NA,ntags)
	}
	
	method <- match.arg(method,c("auto","bin.spline","bin.loess","power","spline"))
	if(method=="auto"){
		if(ntags < 200) {
			method <- "power"
		} else {
			method <- "bin.spline"
		}
	}
	switch(method,
		bin.spline=dispBinTrend(y, design, offset=offset, method.trend="spline", ...),
		bin.loess=dispBinTrend(y, design, offset=offset, method.trend="loess", ...),
		power=dispCoxReidPowerTrend(y, design, offset=offset, ...),
		spline=dispCoxReidSplineTrend(y, design, offset=offset, ...)
	)
}

