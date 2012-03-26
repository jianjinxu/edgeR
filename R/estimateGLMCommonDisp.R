#  Last modified 26 March 2012

estimateGLMCommonDisp <- function(y, design=NULL, offset=NULL, method="CoxReid", ...) 
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design=NULL, offset=NULL, method="CoxReid", ...)
{
    if( is.null(offset) )
        offset <- getOffset(y)
	d <- estimateGLMCommonDisp(y=y$counts, design=design, offset=offset, method=method, ...)
	y$common.dispersion <- d
	y
}

estimateGLMCommonDisp.default <- function(y, design=NULL, offset=NULL, method="CoxReid", ...)
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
		return(NA)
	}
	method <- match.arg(method, c("CoxReid","Pearson","deviance"))
	switch(method,
		CoxReid=dispCoxReid(y, design=design, offset=offset, ...),
		Pearson=dispPearson(y, design=design, offset=offset, ...),
		deviance=dispDeviance(y, design=design, offset=offset, ...)
	)
}

