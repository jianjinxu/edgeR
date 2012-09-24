glmQLFTest <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, abundance.trend=TRUE)
##	Quasi-likelihood F-tests for DGE glms.
##	Davis McCarthy and Gordon Smyth.
##	Created 18 Feb 2011. Last modified 30 Aug 2012.
{
#	Call glmLRT to get most of the results that we need for the QL F-test calculations
	out <- glmLRT(glmfit, coef=coef, contrast=contrast)

#	Calculate squeezed residual variances (the quasi-likelihood over-dispersion parameter)
	s2 <- glmfit$deviance / glmfit$df.residual
	if( abundance.trend )
		s2.fit <- squeezeVar(s2, df=glmfit$df.residual, covariate=glmfit$abundance)
	else
		s2.fit <- squeezeVar(s2, df=glmfit$df.residual)

#	Compute the QL F-statistic
	F <- out$table$LR / out$df.test / s2.fit$var.post
	df.total <- s2.fit$df.prior+glmfit$df.residual
	df.total <- min(df.total, length(s2)*glmfit$df.residual)
	
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

	out$s2.fit <- s2.fit
	out$df.prior <- s2.fit$df.prior
	out$df.total <- df.total
	out
}

