\name{binomTest}
\alias{binomTest}

\title{Multiple Exact Binomial Tests}
\description{
Compute p-values for differential expression for each tag between two digital libraries,
conditioning on the total count for each tag.
The counts are assumed to be Poisson distributed, so the counts in each group as a proportion of the whole are assumed to binomial.
}

\usage{
binomTest(y1, y2, n1=sum(y1), n2=sum(y2), p=n1/(n1+n2))
}

\arguments{
\item{y1}{integer vector giving counts in first library.
Non-integer values are rounded to the nearest integer.}
\item{y2}{integer vector giving counts in second library.
Of same length as \code{x}.
Non-integer values are rounded to the nearest integer.}
\item{n1}{total number of tags in first library.
Non-integer values are rounded to the nearest integer. Not required if \code{p} is supplied.}
\item{n2}{total number of tags in second library.
Non-integer values are rounded to the nearest integer. Not required if \code{p} is supplied.}
\item{p}{expected proportion of y1 to the total under the null hypothesis.}
}

\details{
Computes an exact two-sided binomial test for each tag.
For each tag, the counts are treated as binomial, conditioning on the total number of counts for that tag.
The null hypothesis is that the expected counts are in the same proportions as the library sizes, i.e., that the binomial probability for the first library is \code{n1/(n1+n2)}.

There are a number of ways to conduct two-sided tests.
This function computes the test equivalent to the shortest possible confidence interval for a given confidence level.

This function is the same as \code{sage.test} in the statmod package.
}

\value{
Numeric vector of p-values.
}

\references{
\url{http://en.wikipedia.org/wiki/Binomial_test}
}

\author{Gordon Smyth}

\seealso{
\code{\link[statmod:sage.test]{sage.test}}, \code{\link{binom.test}}
}

\examples{
binomTest(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
#  Univariate equivalents:
binom.test(5,5+30,p=10000/(10000+15000))$p.value
binom.test(10,10+50,p=10000/(10000+15000))$p.value
}

\keyword{htest}
