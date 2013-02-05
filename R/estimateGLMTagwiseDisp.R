# Created March 2011. Last modified 4 Oct 2012.

estimateGLMTagwiseDisp <- function(y, design=NULL, offset=NULL, dispersion=NULL, prior.df=10, trend=TRUE, span=NULL, AveLogCPM=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design=NULL, offset=NULL, dispersion=NULL, prior.df=10, trend=!is.null(y$trended.dispersion), span=NULL, AveLogCPM=NULL, ...)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(trend) {
		if(is.null(dispersion)) dispersion <- y$trended.dispersion
		if(is.null(dispersion)) stop("No trended.dispersion found in data object. Run estimateGLMTrendedDisp first.")
		if(is.null(AveLogCPM)) AveLogCPM <- y$AveLogCPM
		if(is.null(AveLogCPM)) y$AveLogCPM <- AveLogCPM <- aveLogCPM(y$counts,offset=offset)
	} else {
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No common.dispersion found in data object. Run estimateGLMCommonDisp first.")
	}
	out <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=offset, dispersion=dispersion, trend=trend, prior.df=prior.df, AveLogCPM=AveLogCPM, ...)
	y$prior.df <- prior.df
	y$span <- out$span
	y$tagwise.dispersion <- out$tagwise.dispersion
	y
}

estimateGLMTagwiseDisp.default <- function(y, design=NULL, offset=NULL, dispersion, prior.df=10, trend=TRUE, span=NULL, AveLogCPM=NULL, ...)
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(ntags==0) return(numeric(0))
	nlibs <- ncol(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(rep(NA,ntags))
	}

#	Check span
	if(is.null(span)) if(ntags>10) span <- (10/ntags)^0.23 else span <- 1

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset)

#	Call Cox-Reid grid method
	tagwise.dispersion <- dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, prior.df=prior.df, span=span, AveLogCPM=AveLogCPM, ...)

	list(tagwise.dispersion=tagwise.dispersion,span=span)
}
