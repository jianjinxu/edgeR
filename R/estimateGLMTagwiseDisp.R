# Last modified 23 Apr 2012.

estimateGLMTagwiseDisp <- function(y, design=NULL, offset=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design=NULL, offset=NULL, trend=!is.null(y$trended.dispersion), ...)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(trend)  {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) stop("No trended.dispersion found in data object. Run estimateGLMTrendedDisp first.")
		if(is.null(y$abundance)) y$abundance <- mglmOneGroup(y,offset=offset)
	} else {
		dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No common.dispersion found in data object. Run estimateGLMCommonDisp first.")
	}
	y$tagwise.dispersion <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=offset, dispersion=dispersion, trend=trend, abundance=y$abundance, ...)
	y
}

estimateGLMTagwiseDisp.default <- function(y, design=NULL, offset=NULL, dispersion, trend="TRUE", ...)
{
	y <- as.matrix(y)
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(NA,nrow(y))
	}
	dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, ...)
}
