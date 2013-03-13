# Created March 2011. Last modified 13 March 2013.

estimateGLMTagwiseDisp <- function(y, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design=NULL, offset=NULL, dispersion=NULL, prior.df=10, trend=!is.null(y$trended.dispersion), span=NULL, AveLogCPM=NULL, ...)
{
#	If provided as arguments, offset and AveLogCPM over-rule the values stored in y
	if(!is.null(AveLogCPM)) y$AveLogCPM <- AveLogCPM
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	if(!is.null(offset)) y$offset <- expandAsMatrix(offset,dim(y))

#	Find appropriate dispersion
	if(trend) {
		if(is.null(dispersion)) dispersion <- y$trended.dispersion
		if(is.null(dispersion)) stop("No trended.dispersion found in data object. Run estimateGLMTrendedDisp first.")
	} else {
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No common.dispersion found in data object. Run estimateGLMCommonDisp first.")
	}

	d <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=getOffset(y), dispersion=dispersion, trend=trend, prior.df=prior.df, AveLogCPM=y$AveLogCPM, ...)
	y$prior.df <- prior.df
	y$span <- d$span
	y$tagwise.dispersion <- d$tagwise.dispersion
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
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,lib.size=exp(offset))

#	Call Cox-Reid grid method
	tagwise.dispersion <- dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, prior.df=prior.df, span=span, AveLogCPM=AveLogCPM, ...)

	list(tagwise.dispersion=tagwise.dispersion,span=span)
}
