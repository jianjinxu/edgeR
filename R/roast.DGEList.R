roast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean", gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation, var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999)
#	Rotation gene set testing for RNA-Seq data
#	Yunshun Chen, Gordon Smyth
#	Created 19 Dec 2012. Last revised on 4 Feb 2013
{
#	Check design matrix
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("Need at least two groups, or at least two columns for design matrix")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design")

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

	roast(y=z, index=index, design=design, contrast=contrast, set.statistic=set.statistic, gene.weights=gene.weights, array.weights=array.weights, weights=weights, block=block, correlation=correlation, var.prior=var.prior, df.prior=df.prior, trend.var=trend.var, nrot=nrot)
}



mroast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean", gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation, var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999, adjust.method="BH", midp=TRUE, sort="directional")
#	Rotation gene set testing for RNA-Seq data with multiple sets
#	Yunshun Chen, Gordon Smyth
#	Created 8 Jan 2013
{
#	Check design matrix
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("Need at least two groups, or at least two columns for design matrix")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design")

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
	
	mroast(y=z, index=index, design=design, contrast=contrast, set.statistic=set.statistic, gene.weights=gene.weights, array.weights=array.weights, weights=weights, block=block, correlation=correlation, var.prior=var.prior, df.prior=df.prior, trend.var=trend.var, nrot=nrot, adjust.method=adjust.method, midp=midp, sort=sort)
}

