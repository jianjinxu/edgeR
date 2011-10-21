# Last modified 21 Oct 2011.

estimateGLMTagwiseDisp <- function(y, design, offset=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design, offset=NULL, trend=!is.null(y$trended.dispersion), ...)
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

estimateGLMTagwiseDisp.default <- function(y, design, offset=NULL, dispersion, trend="TRUE", ...)
{
	dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, ...)
}
