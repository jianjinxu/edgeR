### Methods for fitting models and testing significance using quasi-likelihood

glmQLFTest <- function(y, glmfit, coef=ncol(glmfit$design), contrast=NULL, abundance.trend=TRUE)
    ##    Tagwise quasi-likelihood F-tests for DGEGLM
    ##    Davis McCarthy and Gordon Smyth.
    ##    Created 18 Feb 2011. Last modified 16 Jan 2012.
{
    y.mat <- as.matrix(y)

    ##    Call glmLRT to get most of the results that we need for the QL F-test calculations
    out.lrt <- glmLRT(y, glmfit, coef=coef, contrast=contrast)

    ##    Calculate squeezed sigma-squared values (the quasi-likelihood parameter)
    s2 <- glmfit$deviance / glmfit$df.residual
    if( abundance.trend )
        s2.fit <- squeezeVar(s2, df=glmfit$df.residual, covariate=glmfit$abundance)
    else
        s2.fit <- squeezeVar(s2, df=glmfit$df.residual)

    ##    Compute the QL F-statistic from the likelihood ratio (LR) statistics, degreees of freedom (df) and QL parameter (sigma-squared)
    LR <- out.lrt$table$LR
    df <- out.lrt$df
    F <- LR / df / s2.fit$var.post
    df.total <- s2.fit$df.prior+glmfit$df.residual
    
    ##    Compute p-values from the QL F-statistic
    F.pvalue <- pf(F, df1=df, df2=df.total, lower.tail = FALSE, log.p = FALSE)

    out.lrt$table$LR <- out.lrt$table$PValue <- NULL
    out.lrt$table$F <- F
    out.lrt$table$PValue <- F.pvalue

    out.lrt$s2.fit <- s2.fit
    out.lrt$df.total <- df.total
    out.lrt
}

