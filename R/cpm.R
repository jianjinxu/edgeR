cpm <- function(x, ...)
UseMethod("cpm")

cpm.DGEList <- function(x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#	Counts per million for a DGEList
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 1 November 2012
{
	lib.size <- x$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*x$samples$norm.factors
	cpm(x$counts,lib.size=lib.size,log=log,prior.count=prior.count)
}

cpm.default <- function(x, lib.size=NULL, log=FALSE, prior.count=0.25, ...)
#	Counts per million for a matrix
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 1 November 2012
{
	x <- as.matrix(x)
	if(is.null(lib.size)) lib.size <- colSums(x)
	if(log) {
		prior.count.scaled <- lib.size/mean(lib.size)*prior.count
		lib.size <- lib.size+2*prior.count.scaled
	}
	lib.size <- 1e-6*lib.size
	if(log)
		log2(t( (t(x)+prior.count.scaled) / lib.size ))
	else
		t(t(x)/lib.size)
}

rpkm <- function(x, gene.length, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
#	Reads per kilobase of gene length per million reads of sequencing
#	Gordon Smyth
#	Created 1 November 2012. Last modified 11 March 2012.
{
	y <- cpm(x=x,normalized.lib.sizes=normalized.lib.sizes,log=log,prior.count=prior.count)
	gene.length.kb <- gene.length/1000
	if(log)
		y-log2(gene.length.kb)
	else
		y/gene.length.kb
}

