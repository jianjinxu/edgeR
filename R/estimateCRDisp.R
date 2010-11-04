estimateCRDisp <- function(y, design=NULL, offset=0, npts=10, min.disp=0, max.disp=2, nselect=200, rowsum.filter=5, tagwise=FALSE, prior.n=10, trend=FALSE, lib.size=NULL, verbose=TRUE)
	## Function to estimate the common dispersion and tagwise disperisons using CoxReid Adjusted Profile-likelihood
	## The function uses cubic spline interpolation in finding the MLEs.
	## Written by Yunshun Chen, August 2010. Last modified by Yunshun Chen, 05 Oct 2010

{
	if( is(y,"DGEList") ) {
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
	if(is.null(design)) {
		if( is(y, "DGEList") ) {
			design <- model.matrix(~y$samples$group)
			if(verbose)
				cat("Design matrix is being formed from the DGEList object, using y$samples$group.\n")
		}
		else
			stop("No design matrix supplied as an argument with matrix of counts.")
	}
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
		cr.com <- rep(0,ntags)
		for(j in 1:ntags) cr.com[j] <- (.maximize.by.interpolation(spline.pts, smoothy[,j]))^4
	} else {	
		cr.com <- (.maximize.by.interpolation(spline.pts, apl.com))^4
	}
	if(cr.com == min.disp || cr.com == max.disp)	
		warning("Common dispersion not within the selected range. Reset the 'min.disp' or the 'max.disp'.")
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

adjustedProfileLik <- function(dispersion, y, design, offset)
## Function to calculate the adjusted profile-likelihood given dispersion, design matrix and response.
## y is simply a table of counts: rows are genes/tags/transcripts, columns are samples/libraries
## offset needs to be a matrix of offsets of the same dimensions as y
## Yunshun Chen, Gordon Smyth
## Created June 2010. Last modified 3 Nov 2010.
{
	if(any(dim(y)!=dim(offset))) stop("offset must be a matrix with the same dimensions as y, the table of counts.")
	tgw.apl <- rep(0,nrow(y))
	start <- matrix(0,nrow(y),ncol(design))
	start[,1] <- glmNBOneGroup(y,offset,dispersion)
	if(dispersion == 0) {
		f <- poisson()
		for(i in 1:nrow(y)) {
			fit <- glm.fit(design,y[i,],offset=offset[i,],family=f,start=start[i,])
			loglik <- -fit$aic/2+fit$qr$rank
			cr <- sum(log(abs(diag(fit$qr$qr)[1:fit$qr$rank])))
			tgw.apl[i] <- loglik - cr
		}
	} else {
		for(i in 1:nrow(y)) {
			fit <- try(.glmnb.fit(design,y[i,],offset=offset[i,],dispersion=dispersion,tol=1e-5,maxit=30,start=start[i,]))
			if(is(fit,"try-error")) {
				tgw.apl[i] <- NA
			} else {
				mu <- fit$fitted
				loglik <- sum(dnbinom(y[i,],size=1/dispersion,mu=mu,log = TRUE))
				R <- chol(crossprod(design,.vecmat(mu/(1+dispersion*mu),design)))
				cr <- sum(log(abs(diag(R))))
				tgw.apl[i] <- loglik - cr
			}
		}
	}
	tgw.apl
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

.glmnb.fit <- function(X,y,dispersion,offset=0,start=NULL,tol=1e-6,maxit=50,trace=FALSE)
#  Fit negative binomial generalized linear model with log link
#  by Levenberg damped Fisher scoring
#  Yunshun Chen and Gordon Smyth
#  2 November 2010.  Last modified 3 November 2010.
{
#  check input
	X <- as.matrix(X)
	n <- nrow(X)
	p <- ncol(X)
	if(p > n) stop("More columns than rows in X")
	y <- as.vector(y)
	if(n != length(y)) stop("length(y) not equal to nrow(X)")
	if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),deviance=numeric(0)))
	if(!(all(is.finite(y)) || all(is.finite(X)))) stop("All values must be finite and non-missing")
	if(any(y < 0)) stop("y must be non-negative")
	maxy <- max(y)
	if(maxy==0) return(list(coefficients=rep(0,p),fitted.values=rep(0,n),deviance=NA))
	phi <- dispersion

#  starting values
	if(is.null(start)) {
		y1 <- pmax(y,1/6)
		fit <- lm.fit(X,log(y1)-offset)
		beta <- fit$coefficients
		mu <- exp(fit$fitted.values+offset)
	} else {
		beta <- start
		mu <- exp(X %*% beta + offset)
	}

	deviance.nb <- function(y,mu,phi) {
		o <- (y < 1e-14) & (mu < 1e-14)
		if(any(o)) {
			if(all(o)) {
				dev <- 0
			} else {
				y <- y[!o]
				mu <- mu[!o]
			}
		}
		y1 <- pmax(y,1/6)
		2*sum(y*log(y1/mu) + (y+1/phi)*log((mu+1/phi)/(y+1/phi)) )
	}

	dev <- deviance.nb(y,mu,phi)

#	Scoring iteration with Levenberg damping
	iter <- 0
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
	repeat {
		iter <- iter+1

#		information matrix
		v.div.mu <- 1+phi*mu
		XVX <- crossprod(X,.vecmat(mu/v.div.mu,X))
		maxinfo <- max(diag(XVX))
		if(iter==1) {
			lambda <- abs(mean(diag(XVX)))/p
			I <- diag(p)
		}

#		score vector
		dl <- crossprod(X,(y-mu)/v.div.mu)

#		Levenberg damping
		betaold <- beta
		devold <- dev
		lev <- 0
		repeat {
			lev <- lev+1

#			trial step
			R <- chol(XVX + lambda*I)
			dbeta <- backsolve(R,backsolve(R,dl,transpose=TRUE))
			beta <- betaold + dbeta
			mu <- exp(X %*% beta + offset)
			dev <- deviance.nb(y,mu,phi)
			if(dev <= devold || dev/max(mu) < 1e-13) break

#			exit if too much damping
			if(lambda/maxinfo > 1e13) {
				beta <- betaold
				warning("Too much damping - convergence tolerance not achievable")
				break
			}

#			step not successful so increase damping
			lambda <- 2*lambda
			if(trace) cat("Damping increased to",lambda,"\n")
		}

#		iteration output
		if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

#		keep exiting if too much damping
		if(lambda/maxinfo > 1e13) break

#		decrease damping if successful at first try
		if(lev==1) lambda <- lambda/10

#		test for convergence
		if( crossprod(dl,dbeta) < tol || dev/max(mu) < 1e-13) break

#		test for iteration limit
		if(iter > maxit) break
	}

	beta <- drop(beta)
	names(beta) <- colnames(X)
	list(coefficients=beta,fitted.values=as.vector(mu),deviance=dev,iter=iter)
}

.vecmat <- function(v,M) {
#	Multiply the rows of matrix by the elements of a vector,
#	i.e., compute diag(v) %*% M
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[1]) stop("Dimensions do not match")
	v * M
}
