glmQLFTest <- function(y, design=NULL, dispersion=NULL, coef=ncol(glmfit$design), contrast=NULL, abundance.trend=TRUE, robust=FALSE, winsor.tail.p=c(0.05,0.1), plot=FALSE)
#	Quasi-likelihood F-tests for DGE glms.
#	Davis McCarthy and Gordon Smyth.
#	Created 18 Feb 2011. Last modified 21 July 2013.
{
#	Initial fit with trended dispersion
	if(is(y,"DGEList")) {
		if(is.null(dispersion)) {
			dispersion <- y$trended.dispersion
			if(is.null(dispersion)) dispersion <- y$common.dispersion
			if(is.null(dispersion)) dispersion <- 0.05
		}
		glmfit <- glmFit(y,design=design,dispersion=dispersion)
	} else {
		glmfit <- y
		disptype <- attr(glmfit$dispersion,"type")
		if(!is.null(disptype)) if(disptype=="tagwise") stop("glmfit should be computed using trended dispersions, not tagwise")
	}

#	Call glmLRT to get most of the results that we need for the QL F-test calculations
	out <- glmLRT(glmfit, coef=coef, contrast=contrast)
	if(is.null(out$AveLogCPM)) out$AveLogCPM <- aveLogCPM(glmfit$counts)

#	Residual deviances
	df.residual <- glmfit$df.residual

#	Adjust df.residual for fitted values at zero
	zerofit <- (glmfit$fitted.values < 1e-14)
	Q <- qr.Q(qr(glmfit$design))
	h <- rowSums(Q^2)
	dffromzeros <- zerofit %*% (1-h)
	df.residual <- drop(round(df.residual-dffromzeros))

#	Empirical Bayes squeezing of the quasi-likelihood variance factors
	s2 <- glmfit$deviance / df.residual
	s2[df.residual==0] <- 0
	s2 <- pmax(s2,0)
	if(abundance.trend) {
		A <- out$AveLogCPM
	} else {
		A <- NULL
	}
	s2.fit <- squeezeVar(s2,df=df.residual,covariate=A,robust=robust,winsor.tail.p=winsor.tail.p)

#	Plot
	if(plot) {
		if(!abundance.trend) A <- out$AveLogCPM
		plot(A,sqrt(sqrt(s2)),xlab="Average Log2 CPM",ylab="Quarter-Root Mean Deviance",pch=16,cex=0.2)
		o <- order(A)
		points(A[o],sqrt(sqrt(s2.fit$var.post[o])),pch=16,cex=0.2,col="red")
		lines(A[o],sqrt(sqrt(s2.fit$var.prior[o])),col="blue")
		legend("topright",pch=16,col=c("black","red"),legend=c("Raw","Squeezed"))
	}

#	Compute the QL F-statistic
	F <- out$table$LR / out$df.test / s2.fit$var.post
	df.total <- s2.fit$df.prior+df.residual
	max.df.residual <- ncol(glmfit$counts)-ncol(glmfit$design)
	df.total <- pmin(df.total, length(s2)*max.df.residual)

#	Compute p-values from the QL F-statistic
	F.pvalue <- pf(F, df1=out$df.test, df2=df.total, lower.tail=FALSE, log.p=FALSE)

#	Ensure is not more significant than chisquare test
	i <- s2.fit$var.post < 1
	if(any(i)) {
		chisq.pvalue <- pchisq(out$table$LR[i], df=out$df.test[i], lower.tail=FALSE, log.p=FALSE)
		F.pvalue[i] <- pmax(F.pvalue[i],chisq.pvalue)
	}

	out$table$LR <- out$table$PValue <- NULL
	out$table$F <- F
	out$table$PValue <- F.pvalue

	out$df.residual.corrected <- df.residual
	out$s2.fit <- s2.fit
	out$df.prior <- s2.fit$df.prior
	out$df.total <- df.total
	out
}

