exactTest.matrix<-function(y1,y2,mus,r,allZeros=rep(FALSE,nrow(y1))) 
    ## Written by Mark Robinson, last modified by Davis McCarthy, 17 June 2009
    ## A function to calculate P-values using a Fisher-like exact test for the Negative Binomial distribution
    ## y1 and y2 are matrices of counts for two given experimental groups (libraries are assumed to be equal in size - adjusted pseudocounts in the edgeR context)
    ## mus is a vector giving the estimated expected value of the count for each tag under the null hypothesis of no difference between the two groups (i.e. common library size * common concentration)
    ## r is the size parameter for the NB distribution (r = 1/phi) - can be either the same or different for each tag
{
    if(nrow(y1)!=nrow(y2))
        stop("Number of rows of y1 not equal to number of rows of y2\n")
    nrows<-nrow(y1)
    pvals<-rep(1,nrows)
    if(any(allZeros)) {
        pvals[!allZeros] <- Recall(y1[!allZeros,],y2[!allZeros,],mus[!allZeros],r[!allZeros])
    } else {
	v<-cbind(rowSums(y1),rowSums(y2))
	n1<-ncol(y1)
	n2<-ncol(y2)
	if (length(mus)==1) 
            mus<-rep(mus,nrows)
	if (length(r)==1)
            r<-rep(r,nrows)
	N<-ceiling(rowSums(v))
	for (i in 1:length(pvals)) {
            ind<-0:N[i]
            p.top<-dnbinom(ind,size=n1*r[i],mu=n1*mus[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mus[i])
            p.obs<-dnbinom(round(v[i,1]),size=n1*r[i],mu=n1*mus[i]) * dnbinom(round(v[i,2]),size=n2*r[i],mu=n2*mus[i])
            keep<- p.top<=p.obs
            p.bot<-dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mus[i])
            pvals[i]<-sum(p.top[keep]/p.bot)
	}
    }
    pvals
}
