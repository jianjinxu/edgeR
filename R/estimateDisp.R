##############################################################
########### Weighted Likelihood Empirical Bayes ##############
##############################################################

estimateDisp <- function(y, design=NULL, offset=NULL, prior.df=NULL, trend.method="locfit", span=NULL, grid.length=21, grid.range=c(-10,10), robust=FALSE, winsor.tail.p=c(0.05,0.1), tol=1e-06)
#  Estimating dispersion using weighted conditional likelihood empirical Bayes.
#  Use GLM approach if a design matrix is given, and classic approach otherwise.
#  It calculates a matrix of likelihoods for each gene at a set of dispersion grid points, and then calls WLEB() to do the shrinkage.
#  Yunshun Chen, Gordon Smyth. Created July 2012. Last modified 4 Feb 2013.
{
	if( !is(y,"DGEList") ) stop("y must be a DGEList")
	group <- y$samples$group <- as.factor(y$samples$group)

	trend <- match.arg(trend.method, c("none", "loess", "locfit", "movingave"))
	ntags <- nrow(y$counts)
	nlibs <- ncol(y$counts)

	# Spline points
	spline.pts <- seq(from=grid.range[1],to=grid.range[2],length=grid.length)
	spline.disp <- 0.1 * 2^spline.pts
	grid.vals <- spline.disp/(1+spline.disp)
	l0 <- matrix(0, ntags, grid.length)

	if(is.null(offset)) offset <- getOffset(y)
	AveLogCPM <- aveLogCPM(y)
	offset <- expandAsMatrix(offset, dim(y))


	# Classic edgeR
	if(is.null(design)){
		# One group
		if(length(levels(group))==1)
			design <- matrix(1,nlibs,1)
		else
			design <- model.matrix(~group)
		if( all(tabulate(group)<=1) ) {
			warning("There is no replication, setting dispersion to NA.")
			y$common.dispersion <- NA
			return(y)
		}
		pseudo.obj <- y

		q2q.out <- equalizeLibSizes(y, dispersion=0.01)
		pseudo.obj$counts <- q2q.out$pseudo
		ysplit <- splitIntoGroups(pseudo.obj)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=ysplit, der=0)
		delta <- delta$maximum
		disp <- delta/(1-delta)

		q2q.out <- equalizeLibSizes(y,dispersion=disp)
		pseudo.obj$counts <- q2q.out$pseudo
		ysplit <- splitIntoGroups(pseudo.obj)
	
		for(j in 1:grid.length) for(i in 1:length(ysplit)) 
			l0[,j] <- condLogLikDerDelta(ysplit[[i]], grid.vals[j], der=0) + l0[,j]
	}
	# GLM edgeR
	else {
		design <- as.matrix(design)
		if(ncol(design) >= ncol(y$counts)) {
			warning("No residual df: setting dispersion to NA")
			y$common.dispersion <- NA
			return(y)
		}
		for(i in 1:grid.length)
			l0[,i] <- adjustedProfileLik(spline.disp[i], y=y$counts, design=design, offset=offset)

	}

	out.1 <- WLEB(theta=spline.pts, loglik=l0, covariate=AveLogCPM, trend.method=trend.method, span=span, individual=FALSE, m0.out=TRUE)

	y$common.dispersion <- 0.1 * 2^out.1$overall
	y$trended.dispersion <- 0.1 * 2^out.1$trend
	y$trend.method <- trend.method
	y$AveLogCPM <- AveLogCPM
	y$span <- out.1$span

	# Calculate prior.df
	if(is.null(prior.df)){
		glmfit <- glmFit(y, design, dispersion=y$trended.dispersion, prior.count=0)

		# Residual deviances
		df.residual <- glmfit$df.residual

		# Adjust df.residual for fitted values at zero
		zerofit <- (glmfit$fitted.values < 1e-14)
		Q <- qr.Q(qr(glmfit$design))
		h <- rowSums(Q^2)
		dffromzeros <- zerofit %*% (1-h)
		df.residual <- drop(round(df.residual-dffromzeros))

		# Empirical Bayes squeezing of the quasi-likelihood variance factors
		s2 <- glmfit$deviance / df.residual
		s2[df.residual==0] <- 0
		s2 <- pmax(s2,0)
		s2.fit <- squeezeVar(s2, df=df.residual, covariate=AveLogCPM, robust=robust, winsor.tail.p=winsor.tail.p)

		prior.df <- s2.fit$df.prior
	}
	ncoefs <- ncol(design)
	prior.n <- prior.df/(nlibs-ncoefs)
			
	# Protecting against infinite prior.n's; otherwise, interpolation of a matrix of Inf values will give the smallest value.
	if(!robust){
	# scalar prior.n
		if (prior.n > 1e6) { 
			if (trend.method!='none') { 
				y$tagwise.dispersion <- y$trended.dispersion
			} else {
				y$tagwise.dispersion <- rep(y$common.dispersion, ntags)
			}
		} else {
			out.2 <- WLEB(theta=spline.pts, loglik=l0, prior.n=prior.n, covariate=AveLogCPM, 
				trend.method=trend.method, span=span, overall=FALSE, trend=FALSE, m0=out.1$shared.loglik)
			y$tagwise.dispersion <- 0.1 * 2^out.2$individual	
		}

	} else {
	# vector prior.n
		i <- prior.n > 1e6
		y$tagwise.dispersion <- rep(y$common.dispersion, ntags)
		if (trend.method!='none') { 
			y$tagwise.dispersion[i] <- y$trended.dispersion[i]
		} 
		if(sum(!i)!=0){
		# Make sure that there are still some genes with finite prior.df
			out.2 <- WLEB(theta=spline.pts, loglik=l0[!i,], prior.n=prior.n[!i], covariate=AveLogCPM[!i], 
				trend.method=trend.method, span=span, overall=FALSE, trend=FALSE, m0=out.1$shared.loglik[!i,])
			y$tagwise.dispersion[!i] <- 0.1 * 2^out.2$individual	
		}
	}
	y$prior.df <- prior.df
	y$prior.n <- prior.n
	y
}



WLEB <- function(theta, loglik, prior.n=5, covariate=NULL, trend.method="locfit", span=NULL, 
	overall=TRUE, trend=TRUE, individual=TRUE, m0=NULL, m0.out=FALSE)
#  Weighted likelihood empirical Bayes for estimating a parameter vector theta
#  given log-likelihood values on a grid of theta values
#  Yunshun Chen, Gordon Smyth
#	Created July 2012. Last modified 24 October 2012.
{
#	Check loglik
	loglik <- as.matrix(loglik)
	ntheta <- ncol(loglik)
	ntags <- nrow(loglik)

#	Check covariate and trend
	if(is.null(covariate))
		trend.method <- "none"
	else
		trend.method <- match.arg(trend.method, c("none", "loess", "locfit", "movingave"))

#	Set span
	if(is.null(span)) if(ntags<=50) span <- 1 else span <- 0.25+0.75*(50/ntags)^0.5

#	Output	
	out <- list()
	out$span <- span

#	overall prior
	if(overall)
		out$overall <- maximizeInterpolant(theta, matrix(colSums(loglik), nrow=1))

#	trended prior
	if(is.null(m0))
	m0 <- switch(trend.method,
		"movingave" = {
			o <- order(covariate)
			oo <- order(o)
			movingAverageByCol(loglik[o,], width=floor(span*ntags))[oo,]
		},
		"loess" = loessByCol(loglik, covariate, span=span)$fitted.values,
		"locfit" = locfitByCol(loglik, covariate, span=span, degree=0),
		"none" = matrix(colMeans(loglik), ntags, length(theta), byrow=TRUE)
	)

	if(trend)
		out$trend <- maximizeInterpolant(theta, m0)

#	weighted empirical Bayes posterior estimates
	if(individual){
		prior.n <- expandAsMatrix(as.vector(prior.n), dim(m0))
		l0a <- loglik + prior.n*m0
		out$individual <- maximizeInterpolant(theta, l0a)
	}

	if(m0.out) out$shared.loglik <- m0

	out
}
