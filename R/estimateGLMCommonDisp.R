estimateGLMCommonDisp <- function(y, design, method="CoxReid", ...) 
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design, method="CoxReid", ...)
{
	d <- estimateGLMCommonDisp(y=y$counts, design=design, offset=getOffsets(y), method=method, ...)
	y$common.dispersion <- d
	y
}

estimateGLMCommonDisp.default <- function(y, design, method="CoxReid", offset=NULL, ...)
{
	y <- as.matrix(y)
	method <- match.arg(method, c("CoxReid","Pearson","deviance"))
	switch(method,
		CoxReid=dispCoxReid(y, design, offset=offset, ...),
		Pearson=dispPearson(y, design, offset=offset, ...),
		deviance=dispDeviance(y, design, offset=offset, ...)
	)
}

