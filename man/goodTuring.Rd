\name{goodTuring}
\alias{goodTuring}

\title{Good-Turing Frequency Estimation}

\description{
Non-parametric empirical Bayes estimates of the frequencies of observed (and unobserved) species.
}

\usage{goodTuring(x, plot=FALSE)}

\arguments{ 
	\item{x}{numeric vector of non-negative integers, representing the observed frequency of each species.}
	\item{plot}{logical, whether to plot log-probability (i.e., log frequencies of frequencies)versus log-frequency.}
}

\value{
A list with components
	\item{count}{observed frequencies, i.e., the unique positive values of \code{x}}
	\item{estimated.p}{estimated proportion of species given the count}
	\item{P0}{estimated combined proportion of all undetected species}
	\item{n0}{number of zeros found in \code{x}}
}

\details{
Observed counts are assumed to be Poisson.
Using an non-parametric empirical Bayes strategy, the algorithm evaluates the posterior expectation of each species mean given its observed count.
The posterior means are then converted to proportions.
In the empirical Bayes step, the counts are smoothed by assuming a log-linear relationship between frequencies and frequencies of frequencies.
The basics of the algorithm are from Good (1953).
Gale and Sampson (1995) proposed a simplied algorithm with a rule for switching between the observed and smoothed frequencies, and it is Gale and Sampson's simplified algorithm that is implemented here.
The number of zero values in \code{x} are not used in the algorithm, but is returned by this function.

Sampson gives a C code version on his webpage at
\url{http://www.grsampson.net/RGoodTur.html}
which gives identical results to this function.
}

\references{
Gale, WA, and Sampson, G (1995).
Good-Turing frequency estimation without tears.
\emph{Journal of Quantitative Linguistics} 2, 217-237.
}

\author{Gordon Smyth}

\examples{
#  True means of observed species
lambda <- rnbinom(10000,mu=2,size=1/10)
lambda <- lambda[lambda>1]

#  Oberved frequencies
Ntrue <- length(lambda)
x <- rpois(Ntrue, lambda=lambda)
freq <- goodTuring(x, plot=TRUE)
}

\keyword{models}
