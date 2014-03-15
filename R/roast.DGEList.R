roast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), ...)
#	Rotation gene set testing for RNA-Seq data
#	Yunshun Chen, Gordon Smyth
#	Created 19 Dec 2012. Last revised on 28 Feb 2014
{
#	Check dispersion estimates in y
	dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")

#	Make default design matrix from group factor
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("design not supplied and samples all belong to the same group")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("design matrix must have at least two columns")

#	Check contrast
	if(length(contrast) == 1) {
		u <- rep.int(0, nbeta)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != nbeta) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")

#	Construct null hypothesis design matrix
	QR <- qr(contrast)
	design0 <- t(qr.qty(QR, t(design))[-1, , drop=FALSE])

#	Null hypothesis fit
	fit.null <- glmFit(y, design0, prior.count=0)

#	Quantile residuals from null fit
	y <- zscoreNBinom(y$counts, mu=fit.null$fitted.values, size=1/dispersion)

	NextMethod("roast")
}



mroast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), ...)
#	Rotation gene set testing for RNA-Seq data with multiple sets
#	Yunshun Chen, Gordon Smyth
#	Created 8 Jan 2013.  Last revised 28 Feb 2014.
{
#	Check dispersion estimates in y
	dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")

#	Make default design matrix from group factor
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("design not supplied and samples all belong to the same group")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("design matrix must have at least two columns")

#	Check contrast
	if(length(contrast) == 1) {
		u <- rep.int(0, nbeta)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != nbeta) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")

#	Construct null hypothesis design matrix
	QR <- qr(contrast)
	design0 <- t(qr.qty(QR, t(design))[-1, , drop=FALSE])

#	Null hypothesis fit
	fit.null <- glmFit(y, design0, prior.count=0)

#	Quantile residuals from null fit
	y <- zscoreNBinom(y$counts, mu=fit.null$fitted.values, size=1/dispersion)

	NextMethod("mroast")
}

