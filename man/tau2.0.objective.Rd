\name{tau2.0.objective}

\alias{tau2.0.objective}

\title{Objective function for tau2}

\description{Objective function for tau2 which is used in the rule of how much to squeeze the dispersion parameters towards the common value}

\usage{ 
tau2.0.objective(tau2.0, info.g, score.g) 
}

\arguments{ 

\item{tau2.0}{scalar, value for tau2}

\item{info.g}{observed information for each gene}

\item{score.g}{observed score (first derivative of log-likelihood) for each gene}

}

\value{ scalar, value of objective function at tau2.0 }

\author{Mark Robinson}

\examples{
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
x<-list(data=y,group=rep(1:2,each=2),lib.size=rep(1000:1001,each=2))
scores <- condLogLikDerDelta(y, delta=0.5, der = 1, doSum = TRUE) 
qA <- quantileAdjust(x, alpha = 10, null.hypothesis = TRUE)
exp.inf <- approx.expected.info(x, d=0.5, qA)
sigma2.0.est <- optimize(tau2.0.objective, c(0, 500), info.g = exp.inf, score.g = scores)$min
}
\keyword{algebra}
