camera.DGEList <- function(y, index, design=NULL, contrast=ncol(design), weights=NULL, use.ranks=FALSE, allow.neg.cor=TRUE, trend.var=FALSE, sort=TRUE)
#	Rotation gene set testing for RNA-Seq data accounting for inter-gene correlation
#	Yunshun Chen, Gordon Smyth
#	Created 07 Jan 2013. Last modified 4 Feb 2013.
{
#	Check design matrix
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("Samples all belong to the same group")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("design matrix must have at least two columns")

#	Check dispersion estimates
	dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")

#	Check contrast
	if(length(contrast) == 1) {
		u <- rep.int(0, nbeta)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != nbeta) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")

#	Null design matrix
	QR <- qr(contrast)
	design0 <- t(qr.qty(QR, t(design))[-1, , drop=FALSE])

#	Null fit
	fit.null <- glmFit(y, design0, prior.count=0)

	z <- zscoreNBinom(y$counts, mu=fit.null$fitted.values, size=1/dispersion)

	camera(y=z, index=index, design=design, contrast=contrast, weights=weights, use.ranks=use.ranks, allow.neg.cor=allow.neg.cor, trend.var=trend.var, sort=sort)
}
