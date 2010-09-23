plotSmear <- function (object, pair=NULL, de.tags=NULL, xlab="logConc", ylab="logFC", pch=19, cex=.2, smearWidth=.5, panel.first=grid(), smooth.scatter=FALSE, ...) {
    ## User-level function for creating an MA-plot for DGE data.
    ## Created by Mark Robinson. Last modified by Davis McCarthy, 12 July 2010.
  if ( !(class(object) %in% c("DGEList", "de4DGEList")) ) 
    stop("Currently only supports DGEList/de4DGEList objects as the object argument.")
  levs.group <- levels(object$samples$group)
  if(length(levs.group)==1)
      stop("Cannot produce an MA-plot with only one group. The one group defined is: ",levs.group)
    if (is.null(pair))
    pair <- levs.group[1:2]
  if( !all(pair %in% levs.group) )
    stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
  stopifnot(length(pair)==2)

  cols1 <- pair[1]==object$samples$group
  cols2 <- pair[2]==object$samples$group

  lib.size <- object$samples$lib.size*object$samples$norm.factors
  x <- rowMeans(t( t(object$counts[,cols1])/lib.size[cols1]))
  y <- rowMeans(t( t(object$counts[,cols2])/lib.size[cols2]))
  
  ylab <- paste(ylab, ":", paste(pair[c(2,1)], collapse="-"), paste="")
	i <- match(de.tags,rownames(object$counts))
	i <- i[!is.na(i)]
  
	maPlot( x, y, xlab=xlab, ylab=ylab, pch=pch, cex=cex, smearWidth=smearWidth, de.tags=i, panel.first=panel.first, smooth.scatter=smooth.scatter, ...)
}

