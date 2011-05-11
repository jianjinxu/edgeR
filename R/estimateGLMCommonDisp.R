estimateGLMCommonDisp <- function(y, design, offset=NULL, method="CoxReid", ...) 
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design, offset=NULL, method="CoxReid", ...)
{
    if( is.null(offset) )
        offset <- getOffset(y)
	d <- estimateGLMCommonDisp(y=y$counts, design=design, offset=offset, method=method, ...)
	y$common.dispersion <- d
	y
}

estimateGLMCommonDisp.default <- function(y, design, offset=NULL, method="CoxReid", ...)
{
	y <- as.matrix(y)
	method <- match.arg(method, c("CoxReid","Pearson","deviance"))
	switch(method,
		CoxReid=dispCoxReid(y, design, offset=offset, ...),
		Pearson=dispPearson(y, design, offset=offset, ...),
		deviance=dispDeviance(y, design, offset=offset, ...)
	)
}

