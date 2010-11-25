### NB GLM fitting genewise using glm.fit()

mglmSimple <- function(y, design, dispersion=0, offset=0, weights=NULL)
##	Fit negative binomial generalized linear model for each transcript
##	to a series of digital expression libraries,
##	using genewise calls to stats:::glm.fit().
##	Requires the {MASS} package for the negative.binomial() family
##	Lower-level function. Takes a matrix of counts (y.mat)

##	Davis McCarthy and Gordon Smyth
##	Created 17 August 2010. Last modified 20 Nov 2010.
{
    require(MASS)
    y.mat <- as.matrix(y)
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
            out <- glm.fit(X,z,w,offset=offset.this,family=f)
            coefficients[i,] <- out$coefficients
            fitted.values[i,] <- fitted(out)
            dev[i] <- out$deviance
            df.residual[i] <- out$df.residual
        }
    }
   list(coefficients=coefficients, df.residual=df.residual, deviance=dev, design=design, 
        offset=offset, dispersion=dispersion, weights=weights, fitted.values=fitted.values)
}

