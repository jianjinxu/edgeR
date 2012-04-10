### NB GLM fitting genewise using glm.fit()

mglmSimple <- function(y, design, dispersion=0, offset=0, weights=NULL)
##	Fit negative binomial generalized linear model for each transcript
##	to a series of digital expression libraries,
##	using genewise calls to stats:::glm.fit().
##	Requires the {MASS} package for the negative.binomial() family
##	Lower-level function. Takes a matrix of counts (y)

##	Davis McCarthy and Gordon Smyth
##	Created 17 August 2010. Last modified 10 Apr 2012.
{
#	Check arguments
	require(MASS)
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	design <- as.matrix(design)
	offset <- expandAsMatrix(offset,dim(y))
	if(!is.null(weights)) {
		weights <- expandAsMatrix(weights,dim(y))
		weights[weights <= 0] <- NA
		y[!is.finite(weights)] <- NA
	} else {
		weights <- array(1,dim(y))
	}

#	Define objects in which to store various results from the glm fits
	coefficients <- matrix(NA,nrow=ngenes,ncol=ncol(design))
	fitted.values <- matrix(NA,nrow=ngenes,ncol=nlibs)
	colnames(coefficients) <- colnames(design)
	rownames(coefficients) <- rownames(y)
	dimnames(fitted.values) <- dimnames(y)
	df.residual <- rep(0,ngenes)
	dev <- rep(NA,ngenes)
	error <- converged <- rep(FALSE,ngenes)

#  If common dispersion, then set glm family here
	if(length(dispersion)>1) {
		common.family <- FALSE
		if(length(dispersion)!=ngenes) stop("length(dispersion) should agree with nrow(y)")
	} else {
		common.family <- TRUE
		if(dispersion > 1e-10)
			f <- negative.binomial(link="log",theta=1/dispersion)
		else
			f <- poisson(link="log")
	}

#	Fit a glm to each gene sequentially

	for (i in 1:ngenes) {
		if(!common.family) {
			if(dispersion[i] > 1e-10)
				f <- negative.binomial(link="log",theta=1/dispersion[i])
			else
				f <- poisson(link="log")
			f$aic <- function(y,n,mu,wt,dev) NA
		}

		z <- as.vector(y[i,])
		obs <- is.finite(z)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			z <- z[obs]
			w <- as.vector(weights[i,obs])
			out <- tryCatch(glm.fit(X,z,w,offset=offset[i,obs],family=f),error=function(e) e)
			if(class(out)[1]=="simpleError") {
				error[i] <- TRUE
			} else {
				coefficients[i,] <- out$coefficients
				fitted.values[i,] <- fitted(out)
				dev[i] <- out$deviance
				df.residual[i] <- out$df.residual
				converged[i] <- out$converged
			}
		}
	}
   list(coefficients=coefficients, df.residual=df.residual, deviance=dev, design=design, 
		offset=offset, dispersion=dispersion, weights=weights, fitted.values=fitted.values,
		converged=converged, error=error)
}

