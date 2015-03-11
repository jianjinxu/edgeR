#  FIT QUASI-LIKELIHOOD GENERALIZED LINEAR MODELS

glmQLFit <- function(y, ...)
UseMethod("glmQLFit")

glmQLFit.DGEList <- function(y, design=NULL, dispersion=NULL, offset=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), ...)
# 	Written by Yunshun Chen and Aaron Lun
#	Created 05 November 2014.
{
	if(is.null(dispersion)) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) dispersion <- 0.05
	}
	if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)

	fit <- glmQLFit(y=y$counts, design=design, dispersion=dispersion, offset=offset, lib.size=NULL, abundance.trend=abundance.trend, AveLogCPM=y$AveLogCPM, robust=robust, winsor.tail.p=winsor.tail.p, weights=y$weights, ...)
	fit$samples <- y$samples
	fit$genes <- y$genes
#	fit$prior.df <- y$prior.df
	fit$AveLogCPM <- y$AveLogCPM
	new("DGEGLM",fit)
}

glmQLFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, abundance.trend=TRUE, AveLogCPM=NULL, robust=FALSE, winsor.tail.p=c(0.05, 0.1), ...)
# 	Fits a GLM and computes quasi-likelihood dispersions for each gene.
# 	Written by Yunshun Chen and Aaron Lun, based on code by Davis McCarthy and Gordon Smyth
# 	Created 15 September 2014. Last modified 05 November 2014.
{
#	Check y
	y <- as.matrix(y)
	ntag <- nrow(y)
	nlib <- ncol(y)

#	Check dispersion
	if(is.null(dispersion)) stop("No dispersion values provided.")

#	Check offset and lib.size
	if(is.null(offset)) {
		if(is.null(lib.size)) lib.size <- colSums(y)
		offset <- log(lib.size)
	}
	offset <- expandAsMatrix(offset,dim(y))

	glmfit <- glmFit(y,design=design,dispersion=dispersion,offset=offset,lib.size=lib.size,...)

#	Setting up the abundances.
	if(abundance.trend) {
		if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y) 
		glmfit$AveLogCPM <- AveLogCPM
	} else {
		AveLogCPM <- NULL
	}

#	Adjust df.residual for fitted values at zero
	zerofit <- (glmfit$fitted.values < 1e-4) & (glmfit$counts < 1e-4)
	df.residual <- .residDF(zerofit, design)

#	Empirical Bayes squeezing of the quasi-likelihood variance factors
	s2 <- glmfit$deviance / df.residual
	s2[df.residual==0] <- 0
	s2 <- pmax(s2,0)
	s2.fit <- squeezeVar(s2,df=df.residual,covariate=AveLogCPM,robust=robust,winsor.tail.p=winsor.tail.p)

#	Storing results
	glmfit$df.residual.zeros <- df.residual
	glmfit$df.prior <- s2.fit$df.prior
	glmfit$var.post <- s2.fit$var.post
	glmfit$var.prior <- s2.fit$var.prior
	glmfit
}


glmQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL)
#	Quasi-likelihood F-tests for DGE glms.
#	Davis McCarthy and Gordon Smyth.
#	Created 18 Feb 2011. Last modified 15 Sep 2014.
{
    if(!is(glmfit,"DGEGLM")) stop("glmfit must be an DGEGLM object produced by glmQLFit") 
	if(is.null(glmfit$var.post)) stop("need to run glmQLFit before glmQLFTest") 
	out <- glmLRT(glmfit, coef=coef, contrast=contrast)

#	Compute the QL F-statistic
	F.stat <- out$table$LR / out$df.test / glmfit$var.post
	df.total <- glmfit$df.prior + glmfit$df.residual.zeros
	max.df.residual <- ncol(glmfit$counts)-ncol(glmfit$design)
	df.total <- pmin(df.total, nrow(glmfit)*max.df.residual)

#	Compute p-values from the QL F-statistic
	F.pvalue <- pf(F.stat, df1=out$df.test, df2=df.total, lower.tail=FALSE, log.p=FALSE)

#	Ensure is not more significant than chisquare test
	i <- glmfit$var.post < 1
	if(any(i)) {
		chisq.pvalue <- pchisq(out$table$LR[i], df=out$df.test[i], lower.tail=FALSE, log.p=FALSE)
		F.pvalue[i] <- pmax(F.pvalue[i],chisq.pvalue)
	}

	out$table$LR <- out$table$PValue <- NULL
	out$table$F <- F.stat
	out$table$PValue <- F.pvalue
	out$df.total <- df.total

	out
}

plotQLDisp <- function(glmfit, xlab="Average Log2 CPM", ylab="Quarter-Root Mean Deviance", pch=16, cex=0.2, col.shrunk="red", col.trend="blue", col.raw="black", ...)
# 	Plots the result of QL-based shrinkage.
#	written by Aaron Lun, based on code by Davis McCarthy and Gordon Smyth
#	15 September 2014
{
	A <- glmfit$AveLogCPM
	if(is.null(A)) A <- aveLogCPM(glmfit)
	s2 <- glmfit$deviance / glmfit$df.residual.zeros
	if(is.null(glmfit$var.post)) { stop("need to run glmQLFit before plotQLDisp") }

	plot(A, sqrt(sqrt(s2)),xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col.raw, ...)
	points(A,sqrt(sqrt(glmfit$var.post)),pch=16,cex=0.2,col=col.shrunk)
	if (length(glmfit$var.prior)==1L) { 
		abline(h=sqrt(sqrt(glmfit$var.prior)), col=col.trend)
	} else {
		o <- order(A)
		lines(A[o],sqrt(sqrt(glmfit$var.prior[o])),col=col.trend)
	}

	legend("topright",pch=16,col=c(col.raw,col.shrunk,col.trend),legend=c("Raw","Squeezed", "Trend"))
	invisible(NULL)
}

