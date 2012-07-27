deviances.function <- function(dispersion)
#	Deviance function for multiple GLMs
#	Gordon Smyth
#	23 November 2010. Last modified 26 Nov 2010.
{
	i <- dispersion>0
	if(all(i)) {
#		All Negative binomial
		deviances <- function(y,mu,dispersion) {
			logymu <- log(y/mu)
			logymu[y<1e-14] <- 0
			2*rowSums(y*logymu + (y+1/dispersion)*log((mu+1/dispersion)/(y+1/dispersion)))
		}
	} else {
		if(any(i)) {
#			Some Poisson, some negative binomial
			deviances <- function(y,mu,dispersion) {
				i <- dispersion>0
				f0 <- deviances.function(0)
				f1 <- deviances.function(1)
				dev <- dispersion
				dev[!i] <- f0(y[!i,,drop=FALSE],mu[!i,,drop=FALSE],0)
				dev[i] <- f1(y[i,,drop=FALSE],mu[i,,drop=FALSE],dispersion[i])
				dev
			}
		} else {
#			All Poisson
			deviances <- function(y,mu,dispersion) {
				logymu <- log(y/mu)
				logymu[y<1e-14] <- 0
				2*rowSums(y*logymu-(y-mu))
			}
		}
	}
	deviances
}

######################################################
######### Simple Line Search glm (Multiple) ##########
######################################################

mglmLS <- function(y,design,dispersion=0,offset=0,coef.start=NULL,tol=1e-5,maxit=50,trace=FALSE)
#  Fit the same negative binomial generalized linear model with log link
#  to multipe response vectors
#  by approximate Fisher scoring with simple line search
#  Yunshun Chen and Gordon Smyth
#  12 November 2010.  Revised 27 July 2012.
{
#	Check input
	X <- as.matrix(design)
	ncoef <- ncol(X)
	if(any(y<0)) stop("y must be non-negative")
	if(is.vector(y)) y <- matrix(y, nrow=1)
	ntags <- nrow(y)
	nlibs <- ncol(y)
	phi <- dispersion
	if(any(phi<0)) stop("dispersions must be non-negative")
	if(all(phi==0)) {
		ispoisson <- TRUE
	} else {
		if(any(phi==0)) stop("Cannot mix zero and positive dispersions")
		ispoisson <- FALSE
	}
	phi <- rep(phi,length=ntags)
	offset <- expandAsMatrix(offset,dim(y))

#	Define deviance functions
	deviances <- deviances.function(dispersion)

#	Transform to orthonormal design matrix
	qrX <- qr(X)
	X <- qr.Q(qrX)

	beta <- matrix(0,ntags,ncoef)
	rownames(beta) <- rownames(y)
	colnames(beta) <- colnames(X)
	stepsize <- meanw <- 1/rowMeans(y)+phi

#	Non-iterative solution for all zero case
	nypos <- rowSums(y>0)
	if(any(nypos<1)) {
#		yi <- y[nypos<1,,drop=FALSE]
#		logyi <- log(yi)
#		logyi[yi==0] <- -30
		z <- -30-offset[nypos<1,,drop=FALSE]
		beta[nypos<1,] <- z %*% X
	}
	
#	Index tags still iterating
	i <- nypos >= 1
	ls.fail <- rep(FALSE,ntags)

#	Starting values
	if(any(i))
	if(is.null(coef.start)) {
		z <- log(pmax(y[i,,drop=FALSE],1/6))-offset[i,,drop=FALSE]
#		beta[i,] <- t(qr.coef(qrX,t(z)))
		beta[i,] <- z %*% X
	} else {
		beta[i,] <- coef.start[i,,drop=FALSE]
	}
	mu <- exp(beta %*% t(X) + offset)
	dimnames(mu) <- dimnames(y)

#	Approximate Fisher scoring iteration
	iter <- 0
	if(trace) {
		cat("Iter",iter,"\n")
		cat("Scoring for",sum(i),"tag(s)\n")
#		print(summary(beta))
	}
	while(any(i)) {
		iter <- iter + 1

#		Tagwise test for convergence
		yi <- y[i,,drop=FALSE]
		mui <- mu[i,,drop=FALSE]
		phii <- phi[i]
		z <- (yi-mui)/(1+phii*mui)
#		dbeta <- t(qr.coef(qrX,t(z)))
		dbeta <- z %*% X
		derivbig <- rowMeans(abs(dbeta)) > tol
#		cat("derivbig",derivbig,"\n")
		i[i] <- derivbig
		i[ls.fail] <- FALSE
#		cat("i",i,"\n")
	
		if(iter > maxit) break

		if(trace) {
			cat("Iter",iter,"\n")
			cat("Scoring for",sum(i),"tag(s)\n")
#			print(summary(beta))
		}

#		Subset to data not yet converged
		if(any(!derivbig)) {
			yi <- yi[derivbig,,drop=FALSE]
			mui <- mui[derivbig,,drop=FALSE]
			phii <- phii[derivbig]
			dbeta <- dbeta[derivbig,,drop=FALSE]
		}

#		Current deviance, and prepare for line search
		devi <- deviances(yi,mui,phii)
		betai <- beta[i,,drop=FALSE]
		offseti <- offset[i,,drop=FALSE]
		stepsizei <- stepsize[i]

#		Index tags active in line search
		j <- i[i]
#		cat("j",j,"\n")

#		Line search until deviance is decreased
		iter.ls <- 0
		while(any(j)) {
			iter.ls <- iter.ls + 1
			if(iter.ls > 50) {
#				cat("Line search iteration limit exceeded at iteration ",iter,"\n")
				k <- which(i)[j]
				ls.fail[k] <- TRUE
				i[k] <- FALSE
				break
			}
			betaj <- betai[j,,drop=FALSE] + stepsizei[j]*dbeta[j,,drop=FALSE]
			muj <- exp(betaj %*% t(X) + offseti[j,,drop=FALSE])
			devj <- deviances(yi[j,,drop=FALSE],muj,phii[j])
			decr <- devj < devi[j]
			if(any(decr)) {
				k <- which(i)[j][decr]
				beta[k,] <- betaj[decr,]
				mu[k,] <- exp(beta[k,,drop=FALSE]%*%t(X)+offset[k,,drop=FALSE])
				if(iter.ls==1) stepsize[k] <- 1.2*stepsize[k]
#				print(betaj[decr,])
				j[j] <- !decr
#				cat("j",j,"\n")
			}
			if(trace) {
				if(iter.ls==1 & any(j)) cat("Step halving:",sum(j),"tag(s), ")
			}
			stepsizei[j] <- stepsizei[j]/3
		}
		if(trace) 
			if(iter.ls>1) cat(iter.ls-1,"iteration(s)\n")
	}

	R <- qr.R(qrX)
	beta <- t(solve(R,t(beta)))
	list(coefficients=beta,fitted.values=mu,fail=which(ls.fail),not.converged=which(i))
}
