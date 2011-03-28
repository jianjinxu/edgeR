### NB GLM fitting genewise using statmod:::glmnb.fit()

mglmLevenberg <- function(y, design, dispersion=0, offset=0, start=NULL)
#	Fit negative binomial generalized linear model for each transcript
#	to a series of Seq libraries,
#	using tagwise calls to statmod:::glmnb.fit().
#	Lower-level function. Takes a matrix of counts (y).

#	Gordon Smyth
#	Created 3 March 2011.
#	Last edited by Yunshun Chen on 8 March 2011
{
#	Check arguments
	require(statmod)
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	design <- as.matrix(design)
	if(length(dispersion)<ngenes) dispersion <- rep(dispersion,length.out=ngenes)
	offset <- expandAsMatrix(offset,dim(y))

	if(!is.null(start)) start <- as.matrix(start)

#	Define objects in which to store various results from the glm fits
	coefficients <- matrix(NA,nrow=ngenes,ncol=ncol(design))
	fitted.values <- matrix(NA,nrow=ngenes,ncol=nlibs)
	colnames(coefficients) <- colnames(design)
	rownames(coefficients) <- rownames(y)
	dimnames(fitted.values) <- dimnames(y)
	df.residual <- rep(0,ngenes)
	dev <- rep(NA,ngenes)

#	Fit a glm to each gene sequentially
	for (i in 1:ngenes) {
		z <- as.vector(y[i,])
		obs <- is.finite(z)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			z <- z[obs]
			if(!is.null(start)) {
				start.val <- start[i,]
			} else {
				start.val <- NULL
			}
			out <- glmnb.fit(X=X,y=z,dispersion=dispersion[i],offset=offset[i,],start=start.val) 
			coefficients[i,] <- out$coefficients
			fitted.values[i,] <- fitted(out)
			dev[i] <- out$deviance
			#df.residual[i] <- out$df.residual
		}
	}
   list(coefficients=coefficients, deviance=dev, design=design, 
		offset=offset, dispersion=dispersion, weights=weights, fitted.values=fitted.values)
}

