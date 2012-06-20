estimateSmoothing<-function(object,verbose=TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, March 2010
# A function to estimate a value for alpha, the weight for the weighted conditional log-likelihood, using approximate empirical Bayes methods
# Smoothing parameter is reported as "prior n", that is the weight given to the common conditional log-likelihood in terms of the number of observations that weight is equivalent to. E.g. having prior n = 10 means that the common conditional log-likelihood is given weight equivalent to 10 tag observations.
{
	if (!is(object,"DGEList"))
		stop("Currently supports DGEList objects")
	group<-object$samples$group
	levs.group<-levels(group)
	d<-object$common.dispersion/(object$common.dispersion + 1)  # common delta
	scores<-0
	for(i in 1:length(levs.group)) {
		if ( sum(group==levs.group[i]) > 1) {
			scores<-scores+condLogLikDerDelta(object$pseudo.alt[,group==levs.group[i]],d,der=1)
		}
	}
	exp.inf<-approx.expected.info(object,d,object$pseudo.alt)
        if(any(exp.inf < 0))
            stop("Algorithm to estimate prior.n failed. Recommend that the common dispersion model or a pre-specified value for prior.n (e.g. prior.n=10) are used.") 
	names(exp.inf)<-rownames(object$counts)
	sigma2.0.est<-.odls(scores.g=scores,info.g=exp.inf)
	alpha<-1/(sigma2.0.est*nrow(object$counts)*mean(exp.inf))
	prior.n <- alpha*nrow(object$counts)
	prior.n
}



.odls <- function(scores.g, info.g) {
# Written by Davis McCarthy, March 2010. A function using Newton-Raphson's method to find the root of the function giving the overdispersion of the likelihood scores
  tau02 <- 0
  iter <- 0
   if (mean(scores.g^2/info.g) < 1) {
     warning("Estimate of overdispersion of likelihood scores is strictly zero.\n Returning min val close to zero (1e-08).")
     tau02 <- 1e-08
     return(tau02)
    }
  repeat {
    iter <- iter + 1
    Q <- sum(scores.g^2/info.g/(1+info.g*tau02) - 1)
    derQ <- -sum(scores.g^2/(1+info.g*tau02)^2)
    dif <- Q/derQ
    tau02  <- tau02 - dif
    if (abs(dif/(sqrt(tau02)+.01)) < 1e-06)
      break
    if (iter > 50) {
      warning("Iteration limit exceeded")
      break
    }
  }
  tau02
}
