#  SUBSET DATA SETS

assign("[.DGEList",
function(object, i, j, ...) {
#  Subsetting for DGEList objects
#  Davis McCarthy, Gordon Smyth 
#  24 September 2009.  Last modified 9 Oct 2009.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$counts <- object$counts[,j,drop=FALSE]
			object$samples <- object$samples[j,,drop=FALSE]
			object$samples$group <- as.factor(as.character(object$samples$group))
			#object$samples$lib.size <- object$samples$lib.size[j,drop=FALSE]
			object$pseudo.alt <- object$pseudo.alt[,j,drop=FALSE]
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
                        object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
                        object$infos <- object$infos[i,drop=FALSE]
			object$pseudo.alt <- object$pseudo.alt[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
                        object$all.zeros <- object$all.zeros[i,drop=FALSE]
		} else {
			object$counts <- object$counts[i,j,drop=FALSE]
			object$samples <- object$samples[j,,drop=FALSE]
                        object$samples$group <- as.factor(as.character(object$samples$group))
			object$pseudo.alt <- object$pseudo.alt[i,j,drop=FALSE]
			#object$samples$lib.size <- object$samples$lib.size[j,drop=FALSE]
			object$conc$conc.common <- object$conc$conc.common[i,drop=FALSE]
			object$conc$conc.group <- object$conc$conc.group[i,,drop=FALSE]
                        object$tagwise.dispersion <- object$tagwise.dispersion[i,drop=FALSE]
			object$infos <- object$infos[i,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
                        object$all.zeros <- object$all.zeros[i,drop=FALSE]
		}
	}
	object
})

assign("[.TopTags",
function(object, i, j, ...)
#  Subsetting for TopTags objects
#	Gordon Smyth
#  7 October 2009.  Last modified 9 October 2009.
{
	if(!missing(j)) stop("Subsetting columns not allowed for TopTags object. Try subsetting object$table instead.",call.=FALSE)
	if(!missing(i)) object$table <- object$table[i,,drop=FALSE]
	object
})

