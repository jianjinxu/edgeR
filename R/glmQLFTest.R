glmQLFit <- function(y, design=NULL, dispersion=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05, 0.1), ...)
# 	Fits a GLM and computes quasi-likelihood dispersions for each gene.
# 	Written by Aaron Lun, based on code by Davis McCarthy and Gordon Smyth
# 	Created 15 September 2014. Last modified 17 September 2014.
{
#	Initial fit with trended dispersion
	if(!is(y,"DGEList")) { stop("y must be a DGEList") }
	if(is.null(dispersion)) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) dispersion <- y$common.dispersion
		if(is.null(dispersion)) dispersion <- 0.05
	}
	glmfit <- glmFit(y,design=design,dispersion=dispersion, ...)
	
#	Setting up the abundances.
	A <- NULL
	if(abundance.trend) {
		if(is.null(y$AveLogCPM)) {
			A <- aveLogCPM(y) # For consistency with call from estimateDisp.
		} else {
			A <- y$AveLogCPM
		}
		glmfit$AveLogCPM <- A
	}

#	Adjust df.residual for fitted values at zero
	zerofit <- (glmfit$fitted.values < 1e-4) & (glmfit$counts < 1e-4)
	df.residual <- .residDF(zerofit, design)

#	Empirical Bayes squeezing of the quasi-likelihood variance factors
	s2 <- glmfit$deviance / df.residual
	s2[df.residual==0] <- 0
	s2 <- pmax(s2,0)
	s2.fit <- squeezeVar(s2,df=df.residual,covariate=A,robust=robust,winsor.tail.p=winsor.tail.p)

#	Storing results
	glmfit$df.residual.zeros <- df.residual
	glmfit$s2.fit <- s2.fit
	glmfit$df.prior <- s2.fit$df.prior
	glmfit
}

glmQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL)
#	Quasi-likelihood F-tests for DGE glms.
#	Davis McCarthy and Gordon Smyth.
#	Created 18 Feb 2011. Last modified 15 Sep 2014.
{
    if(!is(glmfit,"DGEGLM")) stop("glmfit must be an DGEGLM object produced by glmQLFit") 
	if(is.null(glmfit$s2.fit)) stop("need to run glmQLFit before glmQLFTest") 
	out <- glmLRT(glmfit, coef=coef, contrast=contrast)

#	Compute the QL F-statistic
	F.stat <- out$table$LR / out$df.test / glmfit$s2.fit$var.post
	df.total <- glmfit$s2.fit$df.prior + glmfit$df.residual.zeros
	max.df.residual <- ncol(glmfit$counts)-ncol(glmfit$design)
	df.total <- pmin(df.total, nrow(glmfit)*max.df.residual)

#	Compute p-values from the QL F-statistic
	F.pvalue <- pf(F.stat, df1=out$df.test, df2=df.total, lower.tail=FALSE, log.p=FALSE)

#	Ensure is not more significant than chisquare test
	i <- glmfit$s2.fit$var.post < 1
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

plotQLDisp <- function(glmfit, xlab="Average Log2 CPM", ylab="Quarter-Root Mean Deviance", pch=16, cex=0.2, 
		col.shrunk="red", col.trend="blue", col.raw="black", ...)
# 	Plots the result of QL-based shrinkage.
#	written by Aaron Lun, based on code by Davis McCarthy and Gordon Smyth
#	15 September 2014
{
	A <- glmfit$AveLogCPM
	if(is.null(A)) A <- aveLogCPM(glmfit)
	s2 <- glmfit$deviance / glmfit$df.residual.zeros
	if(is.null(glmfit$s2.fit)) { stop("need to run glmQLFit before plotQLDisp") }

	plot(A, sqrt(sqrt(s2)),xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col.raw, ...)
	points(A,sqrt(sqrt(glmfit$s2.fit$var.post)),pch=16,cex=0.2,col=col.shrunk)
	if (length(glmfit$s2.fit$var.prior)==1L) { 
		abline(h=sqrt(sqrt(glmfit$s2.fit$var.prior)), col=col.trend)
	} else {
		o <- order(A)
		lines(A[o],sqrt(sqrt(glmfit$s2.fit$var.prior[o])),col=col.trend)
	}

	legend("topright",pch=16,col=c(col.raw,col.shrunk,col.trend),legend=c("Raw","Squeezed", "Trend"))
	invisible(NULL)
}

