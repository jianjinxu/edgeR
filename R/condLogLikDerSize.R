condLogLikDerSize <- function(y, r, der=1L)
# Calculate derivatives of the conditional log-likelihood function l_g{r}
# with respect to r=1/phi (phi is the dispersion parameter)
# der is derivative (0th deriv is the function)
# For a single group of replicate libraries, all of the same total size
# Written by Mark Robinson
{
#	Vector interpreted as matrix of one row, i.e., one gene
	if (is.vector(y)) {
		y <- matrix(y,nrow=1)
	} else {
		y <- as.matrix(y)
	}

	t <- rowSums(y,na.rm=TRUE)
	n <- rowSums(!is.na(y))
	g <- dim(y)[1]

	switch(der+1L,
		rowSums(lgamma(y+r)) + lgamma(n*r) - lgamma(t+n*r) - n*lgamma(r),
		rowSums(digamma(y+r)) + n*digamma(n*r) - n*digamma(t+n*r) - n*digamma(r),
		rowSums(trigamma(y+r)) + n^2*trigamma(n*r) - n^2*trigamma(t+n*r) - n*trigamma(r)
	)
}
