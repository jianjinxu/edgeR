maximizeInterpolant <- function(x,z,maxit=10,eps=1e-7,plot=FALSE)
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

