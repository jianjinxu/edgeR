#  SUBSET DATA SETS

assign("[.DGEList",
function(object, i, j, ...) {
#  Subsetting for DGEList objects
#  Davis McCarthy, Gordon Smyth 
#  24 September 2009.  Last modified 7 May 2012.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$counts <- object$counts[,j,drop=FALSE]
			object$samples <- droplevels(object$samples[j,,drop=FALSE])
			object$pseudo.counts <- object$pseudo.counts[,j,drop=FALSE]
			object$offset <- object$offset[,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object$counts))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$counts <- object$counts[i,,drop=FALSE]
			object$conc$conc.common <- object$conc$conc.common[i,drop=FALSE]
			object$conc$conc.group <- object$conc$conc.group[i,,drop=FALSE]
			object$abundance <- object$abundance[i,drop=FALSE]
			object$trended.dispersion <- object$trended.dispersion[i,drop=FALSE]
			object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
			object$infos <- object$infos[i,drop=FALSE]
			object$pseudo.counts <- object$pseudo.counts[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$all.zeros <- object$all.zeros[i,drop=FALSE]
			object$offset <- object$offset[i,,drop=FALSE]
		} else {
			object$counts <- object$counts[i,j,drop=FALSE]
			object$samples <- droplevels(object$samples[j,,drop=FALSE])
			object$pseudo.counts <- object$pseudo.counts[i,j,drop=FALSE]
			object$conc$conc.common <- object$conc$conc.common[i,drop=FALSE]
			object$conc$conc.group <- object$conc$conc.group[i,,drop=FALSE]
			object$trended.dispersion <- object$trended.dispersion[i,drop=FALSE]
			object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
			object$infos <- object$infos[i,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$all.zeros <- object$all.zeros[i,drop=FALSE]
			object$offset <- object$offset[i,,drop=FALSE]
		}
	}
	object
})


assign("[.DGEGLM",
function(object, i, j, ...)
#  Subsetting for DGEGLM objects
#  Davis McCarthy, Gordon Smyth
#  11 May 2011.  Last modified 8 April 2013.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(!missing(j))
		stop("Subsetting columns not allowed for DGEGLM object. Try subsetting elements of DGEGLM object instead.",call.=FALSE)
	if(!missing(i)) {
		object$coefficients <- object$coefficients[i,,drop=FALSE]
		object$df.residual <- object$df.residual[i,drop=FALSE]
		object$deviance <- object$deviance[i,drop=FALSE]
		object$offset <- object$offset[i,,drop=FALSE]
		object$genes <- object$genes[i,,drop=FALSE]
		object$trended.dispersion <- object$trended.dispersion[i,drop=FALSE]
		object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
		if(length(object$dispersion)>1) object$dispersion <- object$dispersion[i,drop=FALSE]
		object$weights <- object$weights[i,,drop=FALSE]
		object$fitted.values <- object$fitted.values[i,,drop=FALSE]
		object$abundance <- object$abundance[i,drop=FALSE]
	}
	object
})


assign("[.DGEExact",
function(object, i, j, ...)
#  Subsetting for DGEExact objects
#  Davis McCarthy, Gordon Smyth
#  6 October 2010.  Last modified 8 April 2013.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(!missing(j))
		stop("Subsetting columns not allowed for DGEExact object. Try subsetting object$table instead.",call.=FALSE)
	if(!missing(i)) {
		object$table <- object$table[i,,drop=FALSE]
		object$genes <- object$genes[i,,drop=FALSE]
	}
	object
})


assign("[.DGELRT",
function(object, i, j, ...)
#  Subsetting for DGELRT objects
#  Davis McCarthy, Gordon Smyth	
#  6 April 2011.  Last modified 8 April 2013.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(!missing(j))
		stop("Subsetting columns not allowed for DGELRT object. Try subsetting object$table instead.",call.=FALSE)
	if(!missing(i)) {
		object$table <- object$table[i,,drop=FALSE]
		object$genes <- object$genes[i,,drop=FALSE]
		object$abundance <- object$abundance[i,drop=FALSE]
		object$trended.dispersion <- object$trended.dispersion[i,drop=FALSE]
		object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
		object$dispersion <- object$dispersion[i,drop=FALSE]
		object$coefficients <- object$coefficients[i,,drop=FALSE]
	}
	object
})


assign("[.TopTags",
function(object, i, j, ...)
#  Subsetting for TopTags objects
#  Gordon Smyth 
#  7 October 2009. Last modified 8 April 2013.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$table <- object$table[,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object$counts))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$table <- object$table[i,,drop=FALSE]
		} else {
			object$table <- object$table[i,j,drop=FALSE]
		}
	}
	object
})
