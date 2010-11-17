######################################################
######### Simple Line Search glm (Multiple) ##########
######################################################

mglmLS <- function(y,design,dispersion=0,start=NULL,offset=0,tol=1e-5,maxit=50,trace=FALSE)

#  Fit negative binomial generalized linear model with log link
#  by approximate Fisher scoring with simple line search
#  Yunshun Chen and Gordon Smyth
#  12 November 2010.  Revised 17 Nov 2010.
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
	offset <- matrix(offset,ntags,nlibs,byrow=TRUE)

#	Define deviance functions
	if(ispoisson) {
		deviances <- function(y,mu,phi) {
			logymu <- log(y/mu)
			logymu[y<1e-14] <- 0
			2*rowSums(y*logymu-(y-mu))
		}
	} else {
		deviances <- function(y,mu,phi) {
			logymu <- log(y/mu)
			logymu[y<1e-14] <- 0
			2*rowSums(y*logymu + (y+1/phi)*log((mu+1/phi)/(y+1/phi)))
		}
	}

#	Transform to orthonormal design matrix
	qrX <- qr(X)
	X <- qr.Q(qrX)

	beta <- matrix(0,ntags,ncoef)
	rownames(beta) <- rownames(y)
	colnames(beta) <- colnames(X)
	stepsize <- meanw <- 1/rowMeans(y)+phi

#	Non-iterative solutions for low count cases
	nypos <- rowSums(y>0)
	if(any(nypos<2)) {
		yi <- y[nypos<2,,drop=FALSE]
		logyi <- log(yi)
		logyi[yi==0] <- -10
		z <- logyi-offset[nypos<2,,drop=FALSE]
		beta[nypos<2,] <- z %*% X # t(qr.coef(qrX,t(z)))
#		print(beta[nypos<2,])
	}
	
#	Index tags still iterating
	i <- nypos >= 2
#	cat("i",i,"\n")
	ls.fail <- !i

#	Starting values
	if(is.null(start)) {
		z <- log(pmax(y[i,,drop=FALSE],1/6))-offset[i,,drop=FALSE]
#		beta[i,] <- t(qr.coef(qrX,t(z)))
		beta[i,] <- z %*% X
	} else {
		beta[i,] <- start[i,,drop=FALSE]
	}
	mu <- exp(beta %*% t(X) + offset)
	dimnames(mu) <- dimnames(y)

#	Approximate Fisher scoring iteration
	iter <- 0
	if(trace) {
		cat("Iter =",iter,"\n")
		print(summary(beta))
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
			cat("Iter =",iter,"\n")
			cat("Still iterating",sum(i),"\n")
			print(summary(beta))
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
			stepsizei[j] <- stepsizei[j]/3
		}
		if(trace) cat(iter.ls,"line search iterations\n")
	}

	R <- qr.R(qrX)
	beta <- t(solve(R,t(beta)))
	list(coefficients=beta,fitted=mu,fail=which(ls.fail),not.converged=which(i))
}
