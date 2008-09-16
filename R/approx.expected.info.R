
approx.expected.info<-function(object,d,qA,robust=FALSE) {
  g<-unique(object$group); k1<-object$group==g[1]; k2<-object$group==g[2]
  obs.inf<-condLogLikDerDelta(qA$pseudo[,k1],d,der=2,doSum=FALSE)*(-1)+condLogLikDerDelta(qA$pseudo[,k2],d,der=2,doSum=FALSE)*(-1)
  t<-rowSums(qA$pseudo)
  if (robust) {
    require(MASS)
    inf.lm<-rlm(obs.inf~-1+t)
  } else {
    inf.lm<-lm(obs.inf~-1+t)
  }
  exp.inf.approx<-fitted(inf.lm)
  exp.inf.approx
}

