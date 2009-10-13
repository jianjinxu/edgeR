
interpolateHelper<-function(mu,p,r,count.max,verbose=TRUE)
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to produce the pseudo-counts required
{
	if(length(r)==1) { 
		r<-matrix(r,nrow=nrow(mu),ncol=ncol(mu)) 
	}
	if(length(r)==nrow(mu)) {
		r<-outer(r,rep(1,ncol(mu))) 
	}
	if(sum(is.infinite(r)) != 0) { 
		mx<-max(qpois(p,lambda=mu))+1 
	}
	else { 
		mx<-max(qnbinom(p,size=r,mu=mu))+1 
	}
	pseudo<-matrix(,nrow=nrow(p),ncol=ncol(p))
	for(i in 1:nrow(pseudo)) {
		if(sum(is.infinite(r)) != 0) { # If phi=0 for any tag (i.e. r=1/phi is infinite) then use a Poisson model to adjust the counts
			mx<-max(qpois(p[i,],lambda=mu[i,]))+1
			if(is.infinite(mx)) {
				mx<-count.max[i]*1.5
			}
			for(j in 1:ncol(pseudo)) {
				a<-ppois(-1:(mx-1),lambda=mu[i,j])
				pseudo[i,j]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
			}
		}
		else { # If phi is not equal to zero, then use the negative binomial model
			mx<-max(qnbinom(p[i,],size=r[i,],mu=mu[i,]))+1
			if(is.infinite(mx)) {
				mx<-count.max[i]*1.5
			}
			for(j in 1:ncol(pseudo)) {
				a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
				pseudo[i,j]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
			}
		}
		if(verbose) {
			if(i%%1000==0) cat(".") 
		}
	}
	# May not be necessary to use the code below, if code above definitely will not produce NA values
	v<-which(is.na(pseudo),arr.ind=TRUE)
	if(sum(is.infinite(r)) != 0) { 
		mx<-max(qpois(p,lambda=mu))+1 
	}
	else { 
		mx<-max(qnbinom(p,size=r,mu=mu))+1 
	}
	if(is.infinite(mx))
		mx<-max(count.max)*1.5
	if( nrow(v) > 0 ) {
		for(i in 1:nrow(v)) {
			if(sum(is.infinite(r)) != 0) { # Using Poisson distribution if phi=0 for any tag
				for(j in 1:ncol(v)) {
					a<-ppois(-1:(mx-1),lambda=mu[i,j])
					pseudo[v[i,1],v[i,2]]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
				}
			}
			else { # Otherwise use negative binomial
				for(j in 1:ncol(v)) {
					a<-pnbinom(-1:(mx-1),size=r[i,j],mu=mu[i,j])
					pseudo[v[i,1],v[i,2]]<-suppressWarnings(approx(c(0,a+c(diff(a)/2,0)),c(-.5,0:mx),xout=p[i,j])$y)
				}
			}
			if(verbose) cat("-")
		}
	}
	if(verbose) cat("\n")
	return(pseudo)
}


# Original version of the function
.interpolateHelper<-function (mu, p, r, d, verbose = TRUE) 
{
    if (length(r) == 1) {
        r <- matrix(r, nrow = nrow(mu), ncol = ncol(mu))
    }
    if (length(r) == nrow(mu)) {
        r <- outer(r, rep(1, ncol(mu)))
    }
    mx <- max(qnbinom(p, size = r, mu = mu)) + 1
    psudo <- matrix(, nrow = nrow(p), ncol = ncol(p))
    for (i in 1:nrow(psudo)) {
        mx <- max(qnbinom(p[i, ], size = r[i, ], mu = mu[i, ])) + 1
        if (is.infinite(mx)) 
            mx <- max(d[i, ]) * 1.5
        for (j in 1:ncol(psudo)) {
            a <- pnbinom(-1:(mx - 1), size = r[i, j], mu = mu[i,j])
            psudo[i, j] <- suppressWarnings(approx(c(0, a + c(diff(a)/2,0)), c(-0.5, 0:mx), xout = p[i, j])$y)
        }
        if (verbose) {
            if (i%%1000 == 0) 
                cat(".")
        }
    }
    # May not be necessary to use the code below, if code above definitely will not produce NA values
    #v <- which(is.na(psudo), arr.ind = TRUE)
    #mx <- max(qnbinom(p, size = r, mu = mu)) + 1
    #if (is.infinite(mx)) 
    #    mx <- max(d) * 1.5
    #if (nrow(v) > 0) {
    #    for (i in 1:nrow(v)) {
    #        a <- pnbinom(-1:(mx - 1), size = r[i, j], mu = mu[i,j])
    #        psudo[v[i, 1], v[i, 2]] <- suppressWarnings(approx(c(0,a + c(diff(a)/2, 0)), c(-0.5, 0:mx), xout = p[i,j])$y)
    #        if (verbose) 
    #            cat("-")
    #    }
    #}
    if (verbose) 
        cat("\n")
    return(psudo)
}
