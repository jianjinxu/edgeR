plotSmear <- function (object, pair=NULL, de.tags=NULL, xlab="logConc", ylab="logFC", pch=19, cex=.2, smearWidth=.5, panel.first=grid(), smoothScatter=FALSE, ...) {
    ## User-level function for creating an MA-plot for DGE data.
    ## Created by Mark Robinson. Last modified by Davis McCarthy, 12 July 2010.
  if ( !(class(object) %in% c("DGEList", "de4DGEList")) ) 
    stop("Currently only supports DGEList/de4DGEList objects as the object argument.")
  levs.group <- levels(object$samples$group)
  if(length(levs.group)==1)
      stop("Cannot produce an MA-plot with only one group. The one group defined is: ",levs.group)
  if(dim(object$conc$conc.group)[2] < 2)
      stop("Cannot produce an MA-plot with only one group. Dimensions of object$conc$conc.group are ",dim(object$conc$conc.group)[1]," ",dim(object$conc$conc.group)[2],"\n Try defining more groups and executing estimateCommonDisp() before using plotSmear.\n")
  if (is.null(pair))
    pair <- levs.group[1:2]
  if( !all(pair %in% levs.group) )
    stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
  stopifnot(length(pair)==2)

  col1 <- match(pair[1], levs.group)
  col2 <- match(pair[2], levs.group)

  if( !(any(names(object) %in% "conc")) )
    stop("'conc' is not an element of input object.  Try running estimateCommonDisp() first.")
  D <-  object$conc$conc.group

  ylab <- paste(ylab, ":", paste(pair[c(2,1)], collapse="-"), paste="")
	i <- match(de.tags,rownames(object$counts))
	i <- i[!is.na(i)]
  
	maPlot( D[,col1], D[,col2], xlab=xlab, ylab=ylab, pch=pch, cex=cex, smearWidth=smearWidth, de.tags=i, panel.first=panel.first, smoothScatter=smoothScatter, ...)
}
