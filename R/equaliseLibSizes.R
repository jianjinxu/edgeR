equaliseLibSizes <- function(object, phi=0, N=prod(object$lib.size)^(1/ncol(object$data)), null.hypothesis=TRUE, verbose=TRUE)
# Davis McCarthy, July 2009
# A function that simply adjusts the counts for library size for a fixed value of the dispersion parameter
{
	if (length(phi) > 1 & sum(phi==0) > 0)
		stop("Currently cannot support phi with length > 1 containing exact zeroes\n")
	nrows<-nrow(object$data)
	lib.size<-object$lib.size
	group<-object$group
	levs.group<-levels(group)
	y<-splitIntoGroups(object)
	p<-matrix(0,nrow=nrows,ncol=ncol(object$data))
	mu<-matrix(0,nrow=nrows,ncol=ncol(object$data))
	ps<-estimatePs(object,1/phi)
	if (null.hypothesis==TRUE) {
		for(i in 1:length(levs.group)) {
			if(length(phi)==1 & phi==0) { # phi=0 is equivalent to using a Poisson model - here if common dispersion is used and set to phi=0, then a Poisson model is used 
				if(i==1) cat("Using Poisson model to equalise libraries.\n")
				p[,group==levs.group[i]]<-ppois(y[[i]]-1,lambda=outer(ps$p.common,lib.size[group==levs.group[i]]))+dpois(y[[i]],lambda=outer(ps$p.common,lib.size[group==levs.group[i]]))/2
			}
			else { # Otherwise Use negative binomial model
				p[,group==levs.group[i]]<-pnbinom(y[[i]]-1,size=1/phi,mu=outer(ps$p.common,lib.size[group==levs.group[i]]))+dnbinom(y[[i]],size=1/phi,mu=outer(ps$p.common,lib.size[group==levs.group[i]]))/2
			}
			mu[,group==levs.group[i]]<-outer(ps$p.common,rep(N,sum(group==levs.group[i])))
		}
	} else {
		for(i in 1:length(levs.group)) {
			if(length(phi)==1 & phi==0) { # phi=0 is equivalent to using a Poisson model - here if common dispersion is used and set to phi=0, then a Poisson model is used 
				cat("Using Poisson model to equalise libraries.\n")
				p[,group==levs.group[i]]<-ppois(y[[i]]-1,lambda=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))+dpois(y[[i]],lambda=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))/2
			}
			else { # Otherwise use negative binomial model
				p[,group==levs.group[i]]<-pnbinom(y[[i]]-1,size=1/phi,mu=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))+dnbinom(y[[i]],size=1/phi,mu=outer(ps$p.group[,i],lib.size[group==levs.group[i]]))/2
			}
			mu[,group==levs.group[i]]<-outer(ps$p.group[,i],rep(N,sum(group==levs.group[i])))
		}
	}
	count.max<-apply(object$data,1,max)
	pseudo<-interpolateHelper(mu,p,1/phi,count.max,verbose=verbose)
	pseudo[pseudo<0]<-0 
	pseudo
}

qpoispm <- function (p, lambda, lower.tail = TRUE, log.p = FALSE)
#	Poor man's qpois function, avoiding unnecessary Inf values
#	Gordon Smyth
#	31 July 2009
{
	q1 <- qnorm(p, mean=lambda, sd=sqrt(lambda), lower.tail=lower.tail, log.p=log.p)
	q2 <- qgamma(p, shape=lambda, lower.tail=lower.tail, log.p=log.p)
	(q1+q2)/2-0.5
}

qnbinompm <- function (p, mean, dispersion=0, lower.tail = TRUE, log.p = FALSE)
#	Poor man's qnbinom function, avoiding unnecessary Inf values
#	Parametrized in terms of mean and the over-dispersion parameter
#	Gordon Smyth
#	31 July 2009
{
	r <- 1+dispersion*mean
	v <- mean*r
	q1 <- qnorm(p, mean=mean, sd=sqrt(v), lower.tail=lower.tail, log.p=log.p)
	q2 <- qgamma(p, shape=mean/r, scale=r, lower.tail=lower.tail, log.p=log.p)
	(q1+q2)/2-0.5
}

q2qnbinompm <- function (x, input.mean, output.mean, dispersion=0)
#	Poor man's nbinom quantile adjustment function
#	Parametrized in terms of mean and the over-dispersion parameter
#	Gordon Smyth
#	31 July 2009
{
	ri <- 1+dispersion*input.mean
	vi <- input.mean*ri
	ro <- 1+dispersion*output.mean
	vo <- output.mean*ro
	i <- (x >= input.mean)
	j <- !i
	p1 <- p2 <- q1 <- q2 <- x
	if(any(i)) {
		p1[i] <- pnorm(x[i], mean=input.mean[i], sd=sqrt(vi[i]), lower.tail=FALSE, log.p=TRUE)
		p2[i] <- pgamma(x[i], shape=input.mean[i]/ri[i], scale=ri[i], lower.tail=FALSE, log.p=TRUE)
		q1[i] <- qnorm(p1[i], mean=output.mean[i], sd=sqrt(vo[i]), lower.tail=FALSE, log.p=TRUE)
		q2[i] <- qgamma(p2[i], shape=output.mean[i]/ro[i], scale=ro[i], lower.tail=FALSE, log.p=TRUE)
	}
	if(any(j)) {
		p1[j] <- pnorm(x[j], mean=input.mean[j], sd=sqrt(vi[j]), lower.tail=TRUE, log.p=TRUE)
		p2[j] <- pgamma(x[j], shape=input.mean[j]/ri[j], scale=ri[j], lower.tail=TRUE, log.p=TRUE)
		q1[j] <- qnorm(p1[j], mean=output.mean[j], sd=sqrt(vo[j]), lower.tail=TRUE, log.p=TRUE)
		q2[j] <- qgamma(p2[j], shape=output.mean[j]/ro[j], scale=ro[j], lower.tail=TRUE, log.p=TRUE)
	}
	(q1+q2)/2
}
