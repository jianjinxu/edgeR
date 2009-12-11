\name{tau2.0.objective}
\alias{tau2.0.objective}

\title{Objective Function for Tau2}

\description{Objective function for tau2, which is used in the approximate empirical Bayes rule which determines how much to squeeze the dispersion parameters towards the common value. Tau2 is analogous to the prior variance for each tag/gene in an hierarchical model for the estimators of the tag/genewise dispersion parameter, and must be estimated in order to select the smoothing parameter as an approximate EB rule. }

\usage{ 
tau2.0.objective(tau2.0, info.g, score.g) 
}

\arguments{ 

\item{tau2.0}{scalar, value for tau2}

\item{info.g}{numeric vector, observed information for each tag/gene}

\item{score.g}{scalar, observed score (first derivative of log-likelihood) over all tags/genes}

}

\value{ scalar, value of objective function at tau2.0 }

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
x<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(1000:1001,each=2))
scores <- condLogLikDerDelta(y, delta=0.5, der = 1, doSum = TRUE) 
q2q.out<-equalizeLibSizes(x,disp=1,null.hypothesis=TRUE)
exp.inf<-approx.expected.info(x,d=0.5,q2q.out$pseudo)
sigma2.0.est<-optimize(tau2.0.objective,c(0,500),info.g=exp.inf,score.g=scores)$min
}
\keyword{algebra}
