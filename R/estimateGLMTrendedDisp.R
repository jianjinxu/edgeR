estimateGLMTrendedDisp <- function(y, design, offset=NULL, method="bin.loess", ...) 
UseMethod("estimateGLMTrendedDisp")

estimateGLMTrendedDisp.DGEList <- function(y, design, offset=NULL, method="bin.loess", ...)
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

estimateGLMTrendedDisp.default <- function(y, design, offset=NULL, method="bin.loess", ...)
{
	y <- as.matrix(y)
	method <- match.arg(method, c("bin.spline","bin.loess","power","spline"))
	switch(method,
		bin.spline=dispBinTrend(y, design, offset=offset, method.trend="spline", ...),
		bin.loess=dispBinTrend(y, design, offset=offset, method.trend="loess", ...),
		power=dispCoxReidPowerTrend(y, design, offset=offset, ...),
		spline=dispCoxReidSplineTrend(y, design, offset=offset, ...)
	)
}

