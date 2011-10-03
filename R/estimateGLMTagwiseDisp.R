# Last modified 3 Oct 2011.

estimateGLMTagwiseDisp <- function(y, design, offset=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design, offset=NULL, trend=TRUE, ...)
{
	if( is.null(offset) )
		offset <- getOffset(y)

	if(trend==TRUE)  {
		if( is.null(y$trended.dispersion) & is.null(y$common.dispersion) ) stop("method is trend, but DGEList object has a NULL trended.dispersion slot. Run estimateGLMTrendedDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards trended dispersions.\n")
		if( is.null(y$trended.dispersion) & !is.null(y$common.dispersion) ) {
			warning("trend is TRUE, but DGEList object has a NULL trended.dispersion slot. Common dispersion value used for this calculation.\n")
			dispersion <- y$common.dispersion
			trend <- FALSE
		} else {
			dispersion <- y$trended.dispersion
			abundance <- y$abundance
		}
	} else {
		if( is.null(y$common.dispersion) ) stop("DGEList object has a NULL common.dispersion slot. Run estimateGLMCommonDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards a common value.\n")
		dispersion <- rep(y$common.dispersion, nrow(y))
		abundance <- mglmOneGroup(y$counts,offset=offset)
	}
	d <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=offset, dispersion=dispersion, abundance=abundance, trend=trend, ...)
	y$tagwise.dispersion <- d
	y$abundance <- abundance
	y
}

estimateGLMTagwiseDisp.default <- function(y, design, offset=NULL, dispersion, trend="TRUE", ...)
{
	y <- as.matrix(y)
	dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, ...)
}
