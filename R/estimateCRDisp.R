estimateCRDisp <- function(y, design=NULL, offset=NULL, npts=10, rowsum.filter=5, subset=10000, tagwise=FALSE, prior.n=10, min.n=500, lib.size=NULL, disp.trend=NULL, method.trend="binned-spline", verbose=TRUE)
## Estimate NB dispersion by maximizing the CoxReid Adjusted Profile-likelihood
## The function uses cubic spline interpolation in finding the MLEs.
## 'dispCoxReidPowerTrend' or 'dispCoxReidBinTrend' is called to get the trended common dispersions and the boundaries of tagwise dispersions
## Yunshun Chen
## Created August 2010. Last modified 17 Feb 2011.

{
	message("estimateCRDisp is now obsolete, please use one of estimateGLMCommonDisp, estimateGLMTrendedDisp or estimateGLMTagwiseDisp instead")
	if( is(y,"DGEList") ) {
		if(is.null(y$samples$norm.factors))
			y$samples$norm.factors <- rep(1, ncol(y$counts))
		if(is.null(lib.size))
			lib.size <- y$samples$lib.size*y$samples$norm.factors
		y.mat <- y$counts
	} else {
		y.mat <- as.matrix(y)
		if(is.null(lib.size)) {
			lib.size <- colSums(y.mat)
			if(verbose)
				cat("No lib.size supplied, so lib.size is taken as the column sums of the matrix of counts.\n")
		}
	}
	if(is.null(design)) {
		if( is(y, "DGEList") ) {
			design <- model.matrix(~y$samples$group)
			if(verbose)
				cat("Design matrix is being formed from the DGEList object, using y$samples$group.\n")
		}
		else
			stop("No design matrix supplied as an argument with matrix of counts.")
	}
	if(!is.fullrank(design)) stop("design matrix must be full column rank")
	
	ngenes <- nrow(y.mat)
	narrays <- ncol(y.mat)
	tags.used <- rowSums(y.mat) > rowsum.filter
	y.filt <- y.mat[tags.used,]
	ntags <- nrow(y.filt)

	if(is.null(offset)){
	if( is( y, "DGEList") )
		offset.mat <- getOffset(y)
	else
		offset.mat <- log(lib.size)
		offset.mat.filt <- expandAsMatrix(offset.mat, dim(y.filt))
	}else{
		if( length(offset)==length(y.mat) ) {
			offset.mat <- as.matrix(offset, nrow=ngenes, ncol=narrays)
		} else {
			offset.mat <- matrix(0, nrow=ngenes, ncol=narrays)
			if(length(offset)==narrays | length(offset)==1)
				offset.mat <- matrix(offset, nrow=ngenes, ncol=narrays, byrow=TRUE)
			else 
				stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
		}
		offset.mat.filt <- offset.mat[tags.used,]
	}
	abundance <- mglmOneGroup(y.mat, offset=offset.mat)
	
	method.trend <- match.arg(method.trend, c("binned-spline", "binned-loess", "power"))
	if(!is.null(disp.trend)){
		if(length(disp.trend) != ngenes) stop("dispersion trend has to have the same length as the number of genes")
		tagwise <- TRUE
		disp.trend.filt <- disp.trend[tags.used]
		#abundance <- mglmOneGroup(y.filt, offset=offset.mat.filt)
	} else {
		if( method.trend=="binned-spline" )
			trend <- dispBinTrend(y.filt, design, offset=offset.mat.filt, min.n=min.n, method.trend="spline")
		if( method.trend=="binned-loess" )
			trend <- dispBinTrend(y.filt, design, offset=offset.mat.filt, min.n=min.n, method.trend="loess")
		if( method.trend=="power")
			trend <- dispCoxReidPowerTrend(y.filt, design, offset=offset.mat.filt, subset=subset)
		disp.trend.filt <- trend$dispersion
		#abundance <- trend$abundance
	}
	#cat("Trended dispersions = ", disp.trend.filt, "\n")
	#cat("Abundance = ", abundance, "\n")
	#return(list(abundance = abundance, disp = disp.trend.filt))
	spline.pts <- (0:(npts-1))*2/(npts-1) - 1
	spline.disp <- apl.tgw <- smoothy <- matrix(0, nrow=npts, ncol=ntags)
	for(i in 1:npts) spline.disp[i,] <- disp.trend.filt*2^(4*spline.pts[i])

	abundance.rank <- rank(rowSums(y.filt))
	for(i in 1:npts){
		y.apl <- adjustedProfileLik(spline.disp[i,], y.filt, design=design, offset=offset.mat.filt)
		apl.tgw[i,] <- y.apl
		if(tagwise & prior.n!=0){
            fit <- loess(y.apl ~ abundance.rank, span = 0.3, degree = 0, family = "gaussian", iterations = 1)
            smoothy[i,] <- fitted(fit)
		}
	}
	apl.com <- rowSums(apl.tgw)/ntags
	
	cr.com.filt <- disp.trend.filt
	cr.com <- rep(max(cr.com.filt),ngenes)
	cr.com[tags.used] <- cr.com.filt

	if(tagwise){
		cr.tgw.filt <- rep(0, ntags)
		cr.tgw.all <- rep(max(cr.com), ngenes)
		for(j in 1:ntags) cr.tgw.filt[j] <- disp.trend.filt[j]*2^(4*(maximizeInterpolant(spline.pts, apl.tgw[,j]+ prior.n*smoothy[,j])))
		cr.tgw.all[tags.used] <- cr.tgw.filt
	}
	if(is(y,"DGEList")){
		y$design <- design
		y$abundance <- abundance
		y$CR.common.dispersion <- cr.com
		if(tagwise)	y$CR.tagwise.dispersion <- cr.tgw.all
		return(y)
	} else {
        samples <- data.frame(lib.size=lib.size)
        rownames(samples) <- colnames(y.mat)
		if(tagwise){
			new("DGEList",list(samples=y$samples, counts=y$counts, genes=y$genes, design = design, abundance = abundance,
			CR.common.dispersion=cr.com, CR.tagwise.dispersion=cr.tgw.all))
		} else {
			new("DGEList",list(samples=NULL, counts=y.mat, genes=NULL, design = design, abundance = abundance,
			CR.common.dispersion=cr.com))
		}
	}
}

maximizeInterpolant <- function(x,z,maxit=10,eps=1e-7,plot=FALSE)
#	Maximize a function given a table of values
#	by spline interpolation
#	Gordon Smyth
#	26 August 2010. Modified 1 Sept 2010.
{
	n <- length(z)
	imax <- which.max(z)
	r <- range(x)
	x0 <- x[imax]

#	If maximum occurs at end point, return that value
	if(x0==r[1] || x0==r[2]) return(x0)

	f <- splinefun(x,z)
	if(plot) {
		xx <- seq(from=r[1],to=r[2],length=100)
		zz <- f(xx)
		plot(xx,zz,type="l")
		points(x,z)
	}
	x <- x0
	for (iter in 1:maxit) {
		step <- f(x,deriv=1)/f(x,deriv=2)
		x <- x-step
		if(x<r[1] || x>r[2]) {
			warning("Divergence")
			return(x0)
		}
		if(abs(step) < eps) return(x)
	}
	warning("max iterations exceeded")
	x
}

