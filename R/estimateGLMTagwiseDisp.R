estimateGLMTagwiseDisp <- function(y, design, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design, method="trend", ...)
{
    method <- match.arg(method, c("trend","common"))
    if(method=="trend") {
        if( is.null(y$trended.dispersion) ) stop("method is trend, but DGEList object has a NULL trended.dispersion slot. Run estimateGLMTrendedDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards trended dispersions.\n")
        dispersion <- y$trended.dispersion
        abundance <- y$abundance
    }
    if(method=="common") {
        if( is.null(y$trended.dispersion) ) stop("method is common, but DGEList object has a NULL common.dispersion slot. Run estimateGLMCommonDisp on DGEList object before estimateGLMTagwiseDisp to smooth tagwise dispersions towards a common value.\n")
        dispersion <- rep(y$common.dispersion, nrow(y))
        abundance <- mglmOneGroup(y,offset=getOffsets(y))
    }
	d <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=getOffsets(y), dispersion=dispersion, abundance=abundance, ...)
	y$tagwise.dispersion <- d
    y$abundance <- abundance
	y
}

estimateGLMTagwiseDisp.default <- function(y, design, offset, dispersion, ...)
{
	y <- as.matrix(y)
    dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, ...)
}

