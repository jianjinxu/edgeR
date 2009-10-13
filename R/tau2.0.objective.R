tau2.0.objective<-function(tau2.0,info.g,score.g) 
# Written by Mark Robinson
# A function that can be used in an optimization routine to estimate tau0^2	G<-length(info.g)
{
	denom<-info.g*(1+info.g*tau2.0)
	( mean(score.g^2/denom)-1 )^2
}

