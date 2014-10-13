#  DECIDETESTSDGE.R

decideTestsDGE <- function(object,adjust.method="BH",p.value=0.05)
    ##	Accept or reject hypothesis tests across genes and contrasts
    ##	Davis McCarthy
    ##	15 August 2010. Last modified 19 Jan 2012.
{
    if(!is(object,"DGEExact") & !is(object,"DGELRT")) stop("Need DGEExact or DGELRT object") # Expects a DGEExact or DGELRT object
    decideTests(new("MArrayLM", list(p.value=object$table$PValue, coefficients=object$table$logFC)), method="separate", adjust.method=adjust.method, p.value=p.value)
}



