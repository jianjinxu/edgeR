thinCounts <- function(x,prob=0.5)
#	Binomial thinning of counts
#	Gordon Smyth
#	23 March 2011
{
        x[] <- rbinom(length(x),size=x,prob=prob)
        x
}

