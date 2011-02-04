estimateCRDisp <- function(y, design=NULL, offset=0, npts=10, min.disp=0, max.disp=2, nselect=200, rowsum.filter=5, tagwise=FALSE, prior.n=10, trend=FALSE, lib.size=NULL, verbose=TRUE)
## Estimate NB dispersion by maximizing the CoxReid Adjusted Profile-likelihood
## The function uses cubic spline interpolation in finding the MLEs.
## Yunshun Chen
## Created August 2010. Last modified 04 Feb 2011

{
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
				cat("No lib.size supplied, so lib.size is taken as the column sums of the matrix of counts.")
		}
	}
    if( any( abs(y.mat - round(y.mat)) > .Machine$double.eps^0.5 ) )
        stop("Non-integer entries in y (count data matrix) - all entries must integer-valued for estimateCRDisp() to work.") 
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
	if( length(offset)==length(y.mat) ) {
		offset.mat <- as.matrix(offset, nrow=ngenes, ncol=narrays)
	} else {
		offset.mat <- matrix(0, nrow=ngenes, ncol=narrays)
		if(length(offset)==narrays | length(offset)==1)
			offset.mat <- matrix(offset, nrow=ngenes, ncol=narrays, byrow=TRUE)
		else 
			stop("Number of entries in argument 'offset' incompatible with 'y'. Must have length equal to 1 or to the number of entries in the matrix of counts or to the number of columns in the matrix of counts.\n")
	}
	lib.size.mat <- outer(rep(1,nrow(y.mat)), log(lib.size))
	tags.used <- rowSums(y.mat) > rowsum.filter
	y.filt <- y.mat[tags.used,]
	offset.mat.filt <- offset.mat[tags.used,]
	ntags <- nrow(y.filt)
	if(nselect > ntags) {
		nselect <- ntags
	}
	abundance <- rowSums(y.filt)	
	index <- cbind(abundance,c(1:ntags))
	index.order <- index[order(index[,1],decreasing = TRUE),]
	index.select <- floor(seq(1,ntags,ntags/nselect))
	tags.select <- index.order[index.select,][,2]

	offset.mat.select <- offset.mat.filt[tags.select,]
	y.select <- y.filt[tags.select,]
	lib.size.mat.select <- outer(rep(1,nrow(y.select)), log(lib.size))
	lib.size.mat.filt <- outer(rep(1,nrow(y.filt)), log(lib.size))	
	
	lower <- min.disp^(0.25)
	upper <- max.disp^(0.25)
	spline.pts <- lower + (0:(20-1))*(upper-lower)/(20-1)
	spline.disp <- (spline.pts)^4
	apl.com.select <- c()
	apl.tgw.select <- matrix(0, nrow=20, ncol=nselect)
	for(i in 1:20){
		y.apl <- adjustedProfileLik(spline.disp[i], y.select, design=design, offset=offset.mat.select+lib.size.mat.select)
		apl.tgw.select[i,] <- y.apl
	}
	apl.com.select <- rowSums(apl.tgw.select)/nselect
	bound <- max(apl.com.select) - 5/prior.n
	select <- apl.com.select > bound
	min.disp.new <- min(spline.disp[select])
	max.disp.new <- max(spline.disp[select])

	lower.new <- min.disp.new^(0.25)
	upper.new <- max.disp.new^(0.25)
	spline.pts <- lower.new + (0:(npts-1))*(upper.new-lower.new)/(npts-1)
	spline.disp <- (spline.pts)^4
	apl.tgw <- smoothy <- matrix(0, nrow=npts, ncol=ntags)

	abundance.rank <- rank(rowSums(y.filt))
	for(i in 1:npts){
		y.apl <- adjustedProfileLik(spline.disp[i], y.filt, design=design, offset=offset.mat.filt+lib.size.mat.filt)
		apl.tgw[i,] <- y.apl
		if(trend){
			fit <- loess(y.apl ~ abundance.rank, span = 0.3, degree = 0, family = "gaussian", iterations = 1)
			smoothy[i,] <- fitted(fit)
		}
	}
	apl.com <- rowSums(apl.tgw)/ntags
	if(trend){
			cr.com.filt <- rep(0,ntags)
			for(j in 1:ntags) cr.com.filt[j] <- (.maximize.by.interpolation(spline.pts, smoothy[,j]))^4
			cr.com <- rep(max(cr.com.filt),ngenes)
			cr.com[tags.used] <- cr.com.filt
		} else {	
			cr.com <- (.maximize.by.interpolation(spline.pts, apl.com))^4
			if(cr.com == min.disp || cr.com == max.disp)	
			warning("Common dispersion not within the selected range. Reset the 'min.disp' or the 'max.disp'.")
	}
	if(tagwise){
		cr.tgw.filt <- rep(0, ntags)
		if(trend){
			cr.tgw.all <- rep(max(cr.com), ngenes)
			for(j in 1:ntags) cr.tgw.filt[j] <- (.maximize.by.interpolation(spline.pts, apl.tgw[,j]+ prior.n*smoothy[,j]))^4
		} else {
			cr.tgw.all <- rep(cr.com, ngenes)
			for(j in 1:ntags) cr.tgw.filt[j] <- (.maximize.by.interpolation(spline.pts, apl.tgw[,j]+ prior.n*apl.com))^4
		}
		cr.tgw.all[tags.used] <- cr.tgw.filt
	}
	if(is(y,"DGEList")){
		y$design <- design
		y$CR.common.dispersion=cr.com
		if(tagwise)	y$CR.tagwise.dispersion=cr.tgw.all
		return(y)
	} else {
		if(tagwise){
			new("DGEList",list(samples=y$samples, counts=y$counts, genes=y$genes, design = design, 
			CR.common.dispersion=cr.com, CR.tagwise.dispersion=cr.tgw.all))
		} else {
			new("DGEList",list(samples=y$samples, counts=y$counts, genes=y$genes, design = design, 
			CR.common.dispersion=cr.com))
		}
	}
}

.maximize.by.interpolation <- function(x,z,maxit=10,eps=1e-7,plot=FALSE)
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

