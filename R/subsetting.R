#  SUBSET DATA SETS

assign("[.DGEList",
function(object, i, j, keep.lib.sizes=TRUE)
#  Subsetting for DGEList objects
#  24 September 2009.  Last modified 8 Feb 2015.
{  
#	Recognized components
	IJ <- c("counts","pseudo.counts","offset","weights")
	IX <- c("genes")
	JX <- c("samples")
	I  <- c("AveLogCPM","trended.dispersion","tagwise.dispersion","prior.n","prior.df")
#	Obsolete <- c("conc","infos","all.zeros")

	out <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)
	if(!(missing(i) || keep.lib.sizes)) out$samples$lib.size <- colSums(out$counts)
	if(!missing(j)) out$samples$group <- dropEmptyLevels(out$samples$group)
	out
})

assign("[.DGEGLM",
function(object, i, j)
#  Subsetting for DGEGLM objects
#  11 May 2011.  Last modified 09 May 2014.
{
	if(!missing(j)) stop("Subsetting columns not allowed for DGEGLM object.",call.=FALSE)

#	Recognized components
	IJ <- character(0)
	IX <- c("counts","offset","weights","genes","coefficients","fitted.values","unshrunk.coefficients")
	I  <- c("AveLogCPM","trended.dispersion","tagwise.dispersion","prior.n","prior.df","dispersion","df.residual","deviance","iter","failed")
	JX <- c("samples","design")

	subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)
})


assign("[.DGEExact",
function(object, i, j)
#  Subsetting for DGEExact objects
#  Davis McCarthy, Gordon Smyth
#  6 October 2010.  Last modified 11 Dec 2013.
{
	if(!missing(j)) stop("Subsetting columns not allowed for DGEExact objects.",call.=FALSE)
	if(!missing(i)) {
		object$table <- object$table[i,,drop=FALSE]
		object$genes <- object$genes[i,,drop=FALSE]
	}
	object
})


assign("[.DGELRT",
function(object, i, j)
#  Subsetting for DGELRT objects
#  6 April 2011.  Last modified 09 May 2014.
{
	if(!missing(j)) stop("Subsetting columns not allowed for DGELRT object.",call.=FALSE)

#	Recognized components
	IJ <- character(0)
	IX <- c("counts","offset","weights","genes","coefficients","fitted.values","table","unshrunk.coefficients")
	I  <- c("AveLogCPM","trended.dispersion","tagwise.dispersion","prior.n","prior.df","dispersion","df.residual","deviance","iter","failed","df.test")
	JX <- character(0)

	subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)
})


assign("[.TopTags",
function(object, i, j)
#  Subsetting for TopTags objects
#  7 October 2009. Last modified 11 Dec 2013.
{
	if(!missing(i) || !missing(j)) object$table <- object$table[i,j,drop=FALSE]
	object
})
