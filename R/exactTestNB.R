exactTestNB<-function(pseudo,group,pair=1:2,mus,r) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to calculate P-values using a Fisher-like exact test for the Negative Binomial distribution
{
	nrows<-nrow(pseudo)
    group <- as.factor(group)
	levs.group<-levels(group)
	ngroups<-length(levs.group)
	if( sum(pair[1]==levs.group)==0 | sum(pair[2]==levs.group)==0 ) 
		stop("At least one element of given pair is not a group\n")
	y<-splitIntoGroupsPseudo(pseudo,group,pair)
	pvals<-rep(NA,nrows)
	v<-cbind(rowSums(y$y1),rowSums(y$y2))
	n1<-ncol(y$y1)
	n2<-ncol(y$y2)
	if (length(mus)==1) { mus<-rep(mus,nrows) }
	if (length(r)==1) { r<-rep(r,nrows) }
	N<-ceiling(rowSums(v))
	#N<-round(rowSums(y))
	#N<-floor(rowSums(y))
	for (i in 1:length(pvals)) {
		ind<-0:N[i]
		p.top<-dnbinom(ind,size=n1*r[i],mu=n1*mus[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mus[i])
		p.obs<-dnbinom(round(v[i,1]),size=n1*r[i],mu=n1*mus[i]) * dnbinom(round(v[i,2]),size=n2*r[i],mu=n2*mus[i])
		keep<- p.top<=p.obs
		p.bot<-dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mus[i])
		pvals[i]<-sum(p.top[keep]/p.bot)
	} 
	pvals
}