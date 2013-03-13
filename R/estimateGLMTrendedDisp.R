#  Last modified 11 March 2013

estimateGLMTrendedDisp <- function(y, design=NULL, offset=NULL, AveLogCPM=NULL, method="auto", ...) 
UseMethod("estimateGLMTrendedDisp")

estimateGLMTrendedDisp.DGEList <- function(y, design=NULL, offset=NULL, AveLogCPM=NULL, method="auto", ...)
{
#	If provided as arguments, offset and AveLogCPM over-rule the values stored in y
	if(!is.null(AveLogCPM)) y$AveLogCPM <- AveLogCPM
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	if(!is.null(offset)) y$offset <- expandAsMatrix(offset,dim(y))

	d <- estimateGLMTrendedDisp(y=y$counts, design=design, offset=getOffset(y), AveLogCPM=y$AveLogCPM, method=method, ...)
	y$trended.dispersion <- d$dispersion
	y$trend.method <- d$trend.method
	y$bin.dispersion <- d$bin.dispersion
	y$bin.AveLogCPM <- d$bin.AveLogCPM
	y$design <- d$design
	y
}

estimateGLMTrendedDisp.default <- function(y, design=NULL, offset=NULL, AveLogCPM=NULL, method="auto", ...)
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)

#	Check design
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

#	Check offset
	if(is.null(offset)) {
		lib.size <- colSums(y)
		offset <- log(lib.size)
	}
	offset <- expandAsMatrix(offset,dim(y))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,lib.size=exp(offset))

#	Check method	
	method <- match.arg(method,c("auto","bin.spline","bin.loess","power","spline"))
	if(method=="auto"){
		if(ntags < 200) {
			method <- "power"
		} else {
			method <- "bin.spline"
		}
	}

#	Call lower-level function
	trend <- switch(method,
		bin.spline=dispBinTrend(y, design, offset=offset, method.trend="spline", AveLogCPM=AveLogCPM, ...),
		bin.loess=dispBinTrend(y, design, offset=offset, method.trend="loess", AveLogCPM=AveLogCPM, ...),
		power=dispCoxReidPowerTrend(y, design, offset=offset, AveLogCPM=AveLogCPM, ...),
		spline=dispCoxReidSplineTrend(y, design, offset=offset, AveLogCPM=AveLogCPM, ...)
	)

	trend$design <- design
	trend$AveLogCPM <- AveLogCPM
	trend$trend.method <- method
	trend
}

