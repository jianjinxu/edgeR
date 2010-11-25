\name{betaApproxNBTest}
\alias{betaApproxNBTest}

\title{An Approximate Exact Test for Differences between Two Negative Binomial Groups}

\description{Approximate the tail probabilities of a conditional negative binomial exact test of equality of means between groups. }

\usage{ 
betaApproxNBTest(x1, x2, dispersion)
}

\arguments{ 

	    \item{x1}{vector of observed negative binomial variables for group one}

	    \item{x2}{vector of observed negative binomial variables for group two}

	    \item{dispersion}{vector or scalar providing the value of the NB dispersion parameter for each tag to be used for calculating p-values for differences in mean between the two groups.}

}

\value{
Vector of p-values providing the extent of evidence for difference in means between the two groups.

}

\details{
\code{exactTest} is the user-level function for computing p-values for differential expression between groups in DGE data. However, for tags with extremely large counts, the computation of the tail propbabilities of the conditional negative binomial exact test can be unstable. For such tags, the tail probabilities are well approximated by using a transformed beta distribution (Anderson and Boullion, 1972).

}

\references{
Anderson, Dwane E. and Boullion, Thomas L. Homogeneity test for two negative binomial populations. IEEE Transactions on Reliability, Vol. R-21, No. 2, May 1972.
}

\author{Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
x1<-rnbinom(20,size=1,mu=1000)
x2<-rnbinom(20, size=1, mu=1500)
betaApproxNBTest(x1, x2, dispersion=1)
}

\seealso{
Computing p-values for differential expression for each transcript between two (only) digital gene expression libraries can also be done using the \code{sage.test} function in the \code{statmod} package.
}

\keyword{algebra}

