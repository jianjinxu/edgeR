# Written by Mark Robinson, edited by Davis McCarthy, February 2009
# A function to calculate P-values using Fisher's exact test for the Negative Binomial distribution
exactTestNB<-function(pseudo,group,pair,mus,r,verbose=TRUE) {
	nrows<-nrow(pseudo)
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
	if (verbose) cat("Calculating Fisher exact p-values (dot=1000):")
	N<-ceiling(rowSums(pseudo))
	#N<-round(rowSums(y))
	#N<-floor(rowSums(y))
	for (i in 1:length(pvals)) {
		ind<-0:N[i]
		p.top<-dnbinom(ind,size=n1*r[i],mu=n1*mus[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mus[i])
		p.obs<-dnbinom(round(v[i,1]),size=n1*r[i],mu=n1*mus[i]) * dnbinom(round(v[i,2]),size=n2*r[i],mu=n2*mus[i])
		keep<- p.top<=p.obs
		p.bot<-dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mus[i])
		pvals[i]<-sum(p.top[keep]/p.bot)
		if (verbose) { if (i %% 1000 == 0) { cat(".") } }
	} 
	if (verbose) cat("\n")
	pvals
}



# Mark's original version
# A function to calculate P-values using Fisher's exact test for the Negative Binomial distribution
.exactTestNB<-function(y,g,mus,r,verbose=TRUE) {
	gs<-unique(g)
	y1<-y[,g==gs[1]]; if (is.vector(y1)) { y1<-matrix(y1,ncol=1) }; n1<-ncol(y1)
	y2<-y[,g==gs[2]]; if (is.vector(y2)) { y2<-matrix(y2,ncol=1) }; n2<-ncol(y2)
	pvals<-rep(NA,nrow(y1))
	v<-cbind( rowSums(y1), rowSums(y2))
	if (length(mus)==1) { mus<-rep(mus,nrow(y1)) }
	if (length(r)==1) { r<-rep(r,nrow(y1)) }
	if (verbose) cat("Calculating Fisher exact p-values (dot=1000):")
	N<-ceiling(rowSums(y))	#N<-round(rowSums(y))
	#N<-floor(rowSums(y))
	for (i in 1:length(pvals)) {
		ind<-0:N[i]
		p.top<-dnbinom(ind,size=n1*r[i],mu=n1*mus[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mus[i])
		p.obs<-dnbinom(round(v[i,1]),size=n1*r[i],mu=n1*mus[i]) * dnbinom(round(v[i,2]),size=n2*r[i],mu=n2*mus[i])
		keep<- p.top<=p.obs
		p.bot<-dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mus[i])
		pvals[i]<-sum(p.top[keep]/p.bot)
		if (verbose) { if (i %% 1000 == 0) { cat(".") } }
	} 
	if (verbose) cat("\n")
	pvals
}
