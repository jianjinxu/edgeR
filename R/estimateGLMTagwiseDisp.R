# Last modified 23 Sep 2011.

estimateGLMTagwiseDisp <- function(y, design, offset=NULL, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design, offset=NULL, method="trend", ...)
{
	if( is.null(offset) )
		offset <- getOffset(y)
	method <- match.arg(method, c("trend","common"))
	
	if(method=="common") {
		if( is.null(y$common.dispersion) ) stop("method is common, but DGEList object has a NULL common.dispersion slot. Run estimateGLMCommonDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards a common value.\n")
		dispersion <- rep(y$common.dispersion, nrow(y))
		abundance <- mglmOneGroup(y$counts,offset=offset)
	}
	if(method=="trend") {
		if( is.null(y$trended.dispersion) & is.null(y$common.dispersion) ) stop("method is trend, but DGEList object has a NULL trended.dispersion slot. Run estimateGLMTrendedDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards trended dispersions.\n")
		if( is.null(y$trended.dispersion) & !is.null(y$common.dispersion) ) {
			warning("method is trend, but DGEList object has a NULL trended.dispersion slot. Common dispersion value used for this calculation.\n")
			dispersion <- y$common.dispersion
			method <- "common"
		} else {
			dispersion <- y$trended.dispersion
			abundance <- y$abundance
		}
	}	
	
	d <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=offset, dispersion=dispersion, abundance=abundance, method=method, ...)
	y$tagwise.dispersion <- d
	y$abundance <- abundance
	y
}

estimateGLMTagwiseDisp.default <- function(y, design, offset=NULL, dispersion, method="trend", ...)
{
	y <- as.matrix(y)
	dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, method=method, ...)
}
