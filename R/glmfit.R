#  GENERALIZED LINEAR MODELS

glmFit <- function(y, design, dispersion, offset=0, weights=NULL, lib.size=NULL)
    ##	Fit negative binomial generalized linear model for each transcript
    ##  to a series of digital expression libraries
    ##	Davis McCarthy
    ##	Created 17 August 2010. Last modified 7 October 2010.
    ##  Requires the {MASS} package for the negative.binomial() family
    ##  Lower-level function. Takes a matrix of counts (y.mat)
{
    require(MASS)
    if(!is(y,"DGEList")) y <- DGEList(as.matrix(y))

    if(is.null(lib.size)) lib.size <- y$samples$lib.size*y$samples$norm.factors
    y.mat <- y$counts
    samples <- y$samples
    genes <- y$genes

    narrays <- ncol(y.mat)
    design <- as.matrix(design)
    coefficients <- matrix(NA,nrow=nrow(y.mat),ncol=ncol(design))
    fitted.values <- matrix(NA,nrow=nrow(y.mat),ncol=ncol(y.mat))
    colnames(coefficients) <- colnames(design)
    rownames(coefficients) <- rownames(y.mat)
    dimnames(fitted.values) <- dimnames(y.mat)
    if(!is.null(weights)) {
        weights <- asMatrixWeights(weights,dim(y.mat))
        weights[weights <= 0] <- NA
        y.mat[!is.finite(weights)] <- NA
    } else {
        weights <- array(1,dim(y.mat))
    }
    ngenes <- nrow(y.mat)
                                        # Define objects in which to store various results from the glm fits
    df.residual <- rep(0,ngenes)
    dev <- rep(NA,ngenes)
                                            #  Check length of dispersion vector and that dispersion is not equal to zero
    if(length(dispersion)>1)
        common.family <- FALSE
    else {
        common.family <- TRUE
        if(dispersion > 1e-10)
            f <- negative.binomial(link="log",theta=1/dispersion)
        else
            f <- poisson(link="log")
    }
                                        # Fit a glm to each gene sequentially
    if( length(offset)==length(y.mat) ) {
        offset.mat <- as.matrix(offset, nrow=ngenes, ncol=narrays)
        diff.offsets <- TRUE
    } else {
        diff.offsets <- FALSE
        if(length(offset)==narrays)
            offset.vec <- offset
        if(length(offset)==1)
            offset.vec <- rep(offset,narrays)
        else 
            stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
    }
    for (i in 1:ngenes) {
        if(!common.family) {
            if(dispersion[i] > 1e-10)
                f <- negative.binomial(link="log",theta=1/dispersion[i])
            else
                f <- poisson(link="log")
            f$aic <- function(y,n,mu,wt,dev) NA
        }
        z <- as.vector(y.mat[i,])
        obs <- is.finite(z)
        if(diff.offsets)
            offset.this <- offset.mat[i,]
        else
            offset.this <- offset.vec
                                        # Fit the glm to the data for this gene and store certain useful results
        if(sum(obs) > 0) {
            X <- design[obs,,drop=FALSE]
            z <- z[obs]
            w <- as.vector(weights[i,obs])
            offset.this <- offset.this[obs]
            out <- glm.fit(X,z,w,offset=log(lib.size)+offset.this,family=f)
            coefficients[i,] <- out$coefficients
            fitted.values[i,] <- fitted(out)
            dev[i] <- out$deviance
            df.residual[i] <- out$df.residual
        }
    }
   abundance <- log2( rowMeans(t( t(y.mat)/lib.size )) )
   new("DGEGLM",list(coefficients=coefficients, df.residual=df.residual, deviance=dev, design=design, 
                     offset=offset, samples=samples, genes=genes, dispersion=dispersion, 
                     lib.size=lib.size, weights=weights, fitted.values=fitted.values, abundance=abundance))

}




glmLRT <- function(y,glmfit,coef=ncol(glmfit$design))
    ##	Fit negative binomial generalized linear model for each transcript
    ##	Gordon Smyth. Edited by Davis McCarthy.
    ##	Created 1 July 2010. Last modified, 18 August 2010.
{
    require("MASS")
    if(is(y,"DGEList"))
        y.mat <- y$counts 
    else
        y.mat <- as.matrix(y)
    if(!is(glmfit,"DGEGLM"))
        stop("The glmfit argument must be a DGEGLM object for the full model. Run glmFit with the design matrix of the full model before LR testing.\n")
    design <- glmfit$design
    #### Here - define the design0 and design1 matrices ####
    design <- as.matrix(design)
    nbeta <- ncol(design)
    if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
    coef.names <- colnames(design)
    design0 <- design[,-coef,drop=FALSE]
    fit.null <- glmFit(y,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,lib.size=glmfit$lib.size)
    LR <- fit.null$deviance - glmfit$deviance
    LRT.pvalue <- pchisq(LR, df=( fit.null$df.residual - glmfit$df.residual ), lower.tail = FALSE, log.p = FALSE)
    tab <- data.frame(logConc=glmfit$abundance, logFC=glmfit$coefficients[,coef]/log(2), LR.statistic=LR, p.value=LRT.pvalue)
    rownames(tab) <- rownames(y.mat)
    if(is(y,"DGEList")) {
    	y$counts <- NULL
        y$pseudo.alt <- NULL
    	y$table <- tab  #??maybe not
        y$coefficients.full <- glmfit$coefficients
        y$coefficients.null <- fit.null$coefficients
        y$design <- design
        y$dispersion.used <- glmfit$dispersion
    } else {
    	y <- list(table=tab, coefficients.full=glmfit$coefficients, coefficients.null=fit.null$coefficients, design=design, dispersion.used=glmfit$dispersion)
    }
    new("DGELRT",unclass(y))
}



#glmnb.series <- function(y,design0=NULL,design1=NULL,offset=0,weights=NULL,dispersion,fit.full.model=TRUE,lib.size=NULL,scoretest=FALSE)
    ##	Fit negative binomial generalized linear model for each transcript
    ##  to a series of digital expression libraries
    ##	Gordon Smyth
    ##	Created 20 March 2010. Last revised 28 May 2010 by Davis McCarthy.
    ##  Requires the {MASS} package for the negative.binomial() family
    ##  Requires the {statmod} package for glm.scoretest() function
    ## Lower-level function. Takes a matrix of counts (y.mat)
#{
#    require("MASS")
#    require("statmod")
#    if(is(y,"DGEList")) {
#        if(is.null(lib.size))
#            lib.size <- y$samples$lib.size*y$samples$norm.factors
#        y.mat <- y$counts 
#    } else {
#        y.mat <- as.matrix(y)
#        if(is.null(lib.size))
#            lib.size <- colSums(y)
#    }
#    narrays <- ncol(y.mat)
#    if(is.null(design0))
#        design0 <- matrix(1,narrays,1)
#    design0 <- as.matrix(design0)
#    beta0 <- matrix(NA,nrow=nrow(y.mat),ncol=ncol(design0))
#    if(is.null(design1)) {
#        warning("No design matrix supplied for alternative model (design1). Cannot conduct score or LR tests. Only the null model is fitted here.")
#        fit.full.model <- FALSE
#    }
#    else {
#        design1 <- as.matrix(design1)
#        nbeta1 <- ncol(design1)
#    }
#    nbeta <- ncol(design0)
    ## Check to see if there are rows duplicated between design0 and design1
#    for(i in seq_len(nbeta)) {
#        for(j in seq_len(nbeta1)) {
#            if(identical(design0[,i], design1[,j]))
#              stop("One or more columns of design0 are duplicated in design1. The vector/matrix design1 should only consist of the additional column(s) which would be added to design0 in the alternative model, the significance of which we are trying to test.")
#        }
#    }
#    coef.names <- colnames(design0)
#    if(!is.null(weights)) {
#        weights <- asMatrixWeights(weights,dim(y.mat))
#         weights[weights <= 0] <- NA
#        y.mat[!is.finite(weights)] <- NA
#    } else {
#        weights <- array(1,dim(y.mat))
#    }
#    ngenes <- nrow(y.mat)
#    if(fit.full.model) {
#        coef.names1 <- colnames(design1)
#        fc.coefficient <- matrix(NA,ngenes,nbeta1,dimnames=list(rownames(y),coef.names1))
#    } else {
#       fc.coefficient <- NULL
#    }
                                        # Define objects in which to store various results from the glm fits
#    df.residual.null <- df.residual.full <- rep(0,ngenes)
#    scores <- dev.null <- rep(NA,ngenes)
#    if(fit.full.model)
#        dev.full <- df.full <- rep(NA,ngenes)
#    else
#        dev.full <- df.full <- NULL
                                        #  Check length of dispersion vector and that dispersion is not equal to zero
#    if(length(dispersion)>1)
#        common.family <- FALSE
#    else {
#        common.family <- TRUE
#       if(dispersion > 1e-10)
#            f <- negative.binomial(link="log",theta=1/dispersion)
#        else
#            f <- poisson(link="log")
#   }
                                        # Fit a glm to each gene sequentially
#    if( length(offset)==length(y.mat) ) {
#        offset.mat <- as.matrix(offset, nrow=ngenes, ncol=narrays)
#        diff.offsets <- TRUE
#    } else {
#        diff.offsets <- FALSE
#        if(length(offset)==narrays)
#            offset.vec <- offset
#        if(length(offset)==1)
#            offset.vec <- rep(offset,narrays)
#        else 
#            stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
#    }
#    for (i in 1:ngenes) {
#        if(!common.family) {
#            if(dispersion[i] > 1e-10)
#                f <- negative.binomial(link="log",theta=1/dispersion[i])
#            else
#                f <- poisson(link="log")
#            f$aic <- function(y,n,mu,wt,dev) NA
#        }
#       z <- as.vector(y.mat[i,])
#        obs <- is.finite(z)
#        if(diff.offsets)
#            offset.this <- offset.mat[i,]
#        else
#            offset.this <- offset.vec
                                        # Fit the glm to the data for this gene and store certain useful results
#        if(sum(obs) > 0) {
#            X <- design0[obs,,drop=FALSE]
#            z <- z[obs]
#            w <- as.vector(weights[i,obs])
#            offset.this <- offset.this[obs]
#            out <- glm.fit(X,z,w,offset=log(lib.size)+offset.this,family=f)
#            beta0[i,] <- out$coefficients
#           dev.null[i] <- out$deviance
#            df.residual.null[i] <- out$df.residual
#                                        # Fit the full model (if required) - necessary for conducting LRT
#            if(fit.full.model) {
#                if(scoretest)
#                    scores[i] <- glm.scoretest(out, design1[obs,,drop=FALSE], dispersion=1)
#                X.full <- cbind(X,design1[obs,,drop=FALSE])
#                out.full <- glm.fit(X.full,z,w,offset=log(lib.size)+offset.this,family=f)
#                fc.coefficient[i,] <- out.full$coefficients[ncol(X.full)]
#                dev.full[i] <- out.full$deviance
#                df.residual.full[i] <- out.full$df.residual
#            }
#        }
#    }
#    LR <- dev.null - dev.full
#    LRT.pvalue <- pchisq(LR, df=( df.residual.null-df.residual.full ), lower.tail = FALSE, log.p = FALSE)
#    if(scoretest)
#        score.pvalue <-  2*pnorm(abs(scores),lower.tail=FALSE)
#    list(beta0=beta0, score.stats=scores, score.pvalue=score.pvalue, LR.stats=LR, LRT.pvalue=LRT.pvalue, fc.coefficient=fc.coefficient)
#}



