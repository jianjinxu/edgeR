deviances.function <- function(dispersion)
#	Deviance function for multiple GLMs
#	Gordon Smyth
#	23 November 2010. Last modified 22 Nov 2013.
{
	i <- dispersion>0
	if(all(i)) {
#		All Negative binomial
		deviances <- function(y,mu,dispersion,weights=NULL) {
			logymu <- log(y/mu)
			logymu[y<1e-14] <- 0
			if(is.null(weights)) weights <- 1
			2*rowSums((y*logymu + (y+1/dispersion)*log((mu+1/dispersion)/(y+1/dispersion)))*weights)			  
		}
	} else {
		if(any(i)) {
#			Some Poisson, some negative binomial
			deviances <- function(y,mu,dispersion,weights=NULL) {
				i <- dispersion>0
				f0 <- deviances.function(0)
				f1 <- deviances.function(1)
				dev <- dispersion
				dev[!i] <- f0(y[!i,,drop=FALSE],mu[!i,,drop=FALSE],0,weights[!i,,drop=FALSE])
				dev[i] <- f1(y[i,,drop=FALSE],mu[i,,drop=FALSE],dispersion[i],weights[i,,drop=FALSE])
				dev
			}
		} else {
#			All Poisson
			deviances <- function(y,mu,dispersion,weights=NULL) {
				logymu <- log(y/mu)
				logymu[y<1e-14] <- 0
				if(is.null(weights)) weights <- 1
				2*rowSums((y*logymu-(y-mu))*weights)
			}
		}
	}
	deviances
}
