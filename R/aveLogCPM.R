aveLogCPM <- function(y,dispersion=0.05,offset=0)
#	Fit null (single-group) negative-binomial glm with log-link to DGE data
#  Gordon Smyth
#	18 Aug 2010. Last modified 25 July 2012.
{
	y <- as.matrix(y)
	abundance <- mglmOneGroup(y,dispersion=dispersion,offset=offset)
	log1p(exp(abundance+log(1e6)))/log(2)
}
