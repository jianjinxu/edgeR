
tau2.0.objective<-function(tau2.0,info.g,score.g) {
  G<-length(info.g)
  denom<-info.g*(1+info.g*tau2.0)
  ( mean(score.g^2/denom)-1 )^2
}

