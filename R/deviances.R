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
