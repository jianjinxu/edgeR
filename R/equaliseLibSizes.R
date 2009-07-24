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

