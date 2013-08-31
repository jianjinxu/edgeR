\name{edgeRnews}
\title{edgeR News}
\encoding{UTF-8}

\section{Version 3.3.8}{\itemize{
\item
predFC() with design=NULL now uses normalization factors correctly.
However this use of predFC() to compute counts per million is being phased out in favour of cpm().
}}

\section{Version 3.3.5}{\itemize{
\item
Refinement to cutWithMinN() to make the bin numbers more equal in the worst case.

\item
estimateDisp() now creates the design matrix correctly when the design matrix is not given as an argument and there is only one group.  Previously this case gave an error.

\item
Minor edit to glm.h code.
}}

\section{Version 3.3.4}{\itemize{
\item
plotMDS.DGEList now gives a friendly error message when there are fewer than 3 data columns.
}}

\section{Version 3.3.3}{\itemize{
\item
DGEList() accepts NULL as a possible value again for the group, lib.size and norm.factors arguments.
It is treated the same way as a missing argument.
}}

\section{Version 3.3.2}{\itemize{
\item
Update to cutWithMinN() so that it does not fail even when there are many repeated x values.

\item
Refinement to computation for nbins in dispBinTrend.  Now changes more smoothly with the number of genes.  trace argument is retired.

\item
Fixes to calcNormFactors with method="TMM" so that it takes account of lib.size and refCol if these are preset.

\item
Updates to help pages for the data classes.
}}

\section{Version 3.3.1}{\itemize{
\item
Updates to DGEList() and DGEList-class documentation.
Arguments lib.size, group and norm.factors are now set to their defaults in the function definition rather than set to NULL.
}}

\section{Version 3.2.0}{\itemize{

\item
The User's Guide has a new section on between and within subject designs and a new case study on RNA-seq profiling of unrelated Nigerian individuals.
Section 2.9 (item 2) now gives a code example of how to pre-specify the dispersion value.

\item
New functions estimateDisp() and WLEB() to automate the estimation
of common, trended and tagwise dispersions.
The function estimateDisp() provides a simpler alternative pipeline and in principle replaces all the other dispersion estimation functions, for both glms and for classic edgeR.
It can also incorporate automatic estimation of the prior degrees of freedom, and can do this in a robust fashion.

\item
glmLRT() now permits the contrast argument to be a matrix with multiple columns, making the treatment of this argument analogous to that of the coef argument.

\item
glmLRT() now has a new F-test option.
This option takes into account the uncertainty with which the dispersion is estimated and is more conservative than the default chi-square test.

\item
glmQLFTest() has a number of important improvements.
It now has a simpler alternative calling sequence: it can take either a fitted model object as before, or it can take a DGEList object and design matrix and do the model fit itself.
If provided with a fitted model object, it now checks whether the dispersion is of a suitable type (common or trended).
It now optionally produces a plot of the raw and shrunk residual mean deviances versus AveLogCPM.
It now has the option of robustifying the empirical Bayes step.
It now has a more careful calculation of residual df that takes special account of cases where all replicates in a group are identically zero.

\item
The gene set test functions roast(), mroast() and camera() now have methods defined for DGEList data objects.
This facilitates gene set testing and pathway analysis of expression profiles within edgeR.

\item
The default method of plotMDS() for DGEList objects has changed.
The new default forms log-counts-per-million and computes Euclidean distances.
The old method based on BCV-distances is available by setting method="BCV".
The annotation of the plot axes has been improved so that the distance method used is apparent from the plot.

\item
The argument prior.count.total used for shrinking log-fold-changes has been changed to prior.count in various functions throughout the package, and now refers to the average prior.count per observation rather than the total prior count across a transcript.
The treatment of prior.counts has also been changed very slightly in cpm() when log=TRUE.

\item
New function aveLogCPM() to compute the average log count per million for each transcript across all libraries.
This is now used by all functions in the package to set AveLogCPM, which is now the standard measure of abundance. 
The value for AveLogCPM is now computed just once, and not updated when the dispersion is estimated or when a linear model is fitted.
glmFit() now preserves the AveLogCPM vector found in the DGEList object rather than recomputing it.
The use of the old abundance measure is being phased out.

\item
The glm dispersion estimation functions are now much faster.

\item
New function rpkm() to compute reads per kilobase per million (RPKM).

\item
New option method="none" for calcNormFactors().

\item
The default span used by dispBinTrend() has been reduced.

\item
Various improvements to internal C++ code.

\item
Functions binCMLDispersion() and bin.dispersion() have been removed as obsolete.

\item
Bug fix to subsetting for DGEGLM objects.

\item
Bug fix to plotMDS.DGEList to make consistent use of norm.factors.
}}


\section{Version 3.0.0}{\itemize{

\item
New chapter in the User's Guide covering a number of common types of experimental designs, including multiple groups, multiple factors and additive models.
New sections in the User's Guide on clustering and on making tables of read counts.
Many other updates to the User's Guide and to the help pages.

\item
New function edgeRUsersGuide() to open the User's Guide in a pdf viewer.

\item
Many functions have made faster by rewriting the core computations in C++.
This includes adjustedProfileLik(), mglmLevenberg(), maximizeInterpolant() and goodTuring().

\item
New argument verbose for estimateCommonDisp() and estimateGLMCommonDisp().

\item
The trended dispersion methods based on binning and interpolation have been rewritten to give more stable results when the number of genes is not large.

\item
The amount by which the tagwise dispersion estimates are squeezed towards the global value is now specified in estimateTagwiseDisp(), estimateGLMTagwiseDisp() and dispCoxReidInterpolateTagwise() by specifying the prior degrees of freedom prior.df instead of the prior number of samples prior.n.

\item
The weighted likelihood empirical Bayes code has been simplified or developed in a number of ways.
The old functions weightedComLik() and weightedComLikMA() are now removed as no longer required.

\item
The functions estimateSmoothing() and approx.expected.info() have been removed as no longer recommended.

\item
The span used by estimateGLMTagwiseDisp() is now chosen by default as a decreasing function of the number of tags in the dataset.

\item
New method "loess" for the trend argument of estimateTagwiseDisp, with "tricube" now treated as a synonym.

\item
New functions loessByCol() and locfitByCol() for smoothing columns of matrix by non-robust loess curves.
These functions are used in the weighted likelihood empirical Bayes procedures to compute local common likelihood.

\item
glmFit now shrinks the estimated fold-changes towards zero.
The default shrinkage is as for exactTest().

\item
predFC output is now on the natural log scale instead of log2.

\item
mglmLevenberg() is now the default glm fitting algorithm, avoiding the occasional errors that occurred previously with mglmLS().

\item
The arguments of glmLRT() and glmQLFTest() have been simplified so that the argument y, previously the first argument of glmLRT, is no longer required.

\item
glmQLFTest() now ensures that no p-value is smaller than what would be obtained by treating the likelihood ratio test statistic as chisquare.

\item
glmQLFTest() now treats tags with all zero counts in replicate arrays as having zero residual df.

\item
gof() now optionally produces a qq-plot of the genewise goodness of fit statistics.

\item
Argument null.hypothesis removed from equalizeLibSizes().

\item
DGEList no longer outputs a component called all.zeros.

\item
goodTuring() no longer produces a plot.
Instead there is a new function goodTuringPlot() for plotting log-probability versus log-frequency.
goodTuring() has a new argument 'conf' giving the confidence factor for the linear regression approximation.

\item
Added plot.it argument to maPlot().
}}


\section{Version 2.6.0}{\itemize{

\item
edgeR now depends on limma.

\item
Considerable work on the User's Guide.
New case study added on Pathogen inoculated arabidopsis illustrating a two group comparison with batch effects.
All the other case studies have been updated and streamlined.
New section explaining why adjustments for GC content and mappability are not necessary in a differential expression context.

\item
New and more intuitive column headings for topTags() output.
'logFC' is now the first column.
Log-concentration is now replaced by log-counts-per-million ('logCPM').
'PValue'  replaces 'P.Value'.
These column headings are now inserted in the table of results by
exactTest() and glmLRT() instead of being modified by the show method for the TopTags object generated by topTags().
This means that the column names will be correct even when users access the fitted model objects
directly instead of using the show method.

\item
plotSmear() and plotMeanVar() now use logCPM instead of logConc.

\item
New function glmQLFTest() provides quasi-likelihood hypothesis testing using F-tests, as an alternative to likelihood ratio tests using the chisquare distribution.

\item
New functions normalizeChIPtoInput() and calcNormOffsetsforChIP()
for normalization of ChIP-Seq counts relative to input control.

\item
New capabilities for formal shrinkage of the logFC.
exactTest() now incorporates formal shrinkage of the logFC, controlled by argument 'prior.count.total'.
predFC() provides similar shrinkage capability for glms.

\item
estimateCommonDisp() and estimateGLMCommonDisp() now set the dispersion to NA when there is no replication,
instead of setting the dispersion to zero.
This means that users will need to set a dispersion value explicitly to use functions further down the analysis pipeline.

\item
New function estimateTrendedDisp() analogous to estimateGLMTrendedDisp() but for classic edgeR.

\item
The algorithms implemented in estimateTagwiseDisp() now uses fewer grid points but interpolates, similar to estimateGLMTagwiseDisp().

\item
The power trend fitted by dispCoxReidPowerTrend() now includes a positive asymptote.
This greatly improves the fit on real data sets.
This now becomes the default method for estimateGLMTrendedDisp() when the number of genes is less than 200.

\item
New user-friendly function plotBCV() displays estimated dispersions.  

\item
New argument target.size for thinCounts().

\item
New utility functions getDispersion() and zscoreNBinom().

\item
dimnames() methods for DGEExact, DGELRT and TopTags classes.

\item
Function pooledVar() removed as no longer necessary.

\item
Minor fixes to various functions to ensure correct results in special cases.

}}


\section{Version 2.4.0}{\itemize{

\item New function spliceVariants() for detecting alternative exon usage
      from exon-level count data.

\item A choice of rejection regions is now implemented for exactTest(),
      and the default is changed from one based on small probabilities
      to one based on doubling the smaller of the tail probabilities.
      This gives better results than the original conditional test when
      the dispersion is large (especially > 1).  A Beta distribution
      approximation to the tail probability is also implemented when
      the counts are large, making exactTest() much faster and
      less memory hungry.

\item estimateTagwiseDisp() now includes an abundance trend on the
      dispersions by default.
      
\item exactTest() now uses tagwise.dispersion by default if found in the
      object.
    
\item estimateCRDisp() is removed. It is now replaced by
      estimateGLMCommonDisp(), estimateGLMTrendedDisp() and
      estimateGLMTagwiseDisp().
    
\item Changes to glmFit() so that it automatically detects dispersion
      estimates if in data object. It uses tagwise if available, then
      trended, then common.

\item Add getPriorN() to calculate the weight given to the common
      parameter likelihood in order to smooth (or stabilize) the
      dispersion estimates. Used as default for estimateTagwiseDisp and
      estimateGLMTagwiseDisp().
    
\item New function cutWithMinN() used in binning methods.
    
\item glmFit() now S3 generic function, and glmFit() has new method
      argument specifying fitting algorithm.
    
\item DGEGLM objects now subsettable.
    
\item plotMDS.dge() is retired, instead a DGEList method is now defined for 
      plotMDS() in the limma package.  One advantage is that the plot can
      be repeated with different graphical parameters without recomputing
      the distances.  The MDS method is also now much faster.
    
\item Add as.data.frame method for TopTags objects.
    
\item New function cpm() to calculate counts per million conveniently. 
    
\item Adding args to dispCoxReidInterpolateTagwise() to give more access to 
      tuning parameters.
    
\item estimateGLMTagwiseDisp() now uses trended.dispersion by default if 
      trended.dispersion is found.
      
\item Change to glmLRT() to ensure character coefficient argument will work.
    
\item Change to maPlot() so that any really extreme logFCs are brought back 
      to a more reasonable scale.
      
\item estimateGLMCommonDisp() now returns NA when there are no residual
      df rather than returning dispersion of zero.

\item The trend computation of the local common likelihood in
      dispCoxReidInterpolateTagwise() is now based on moving averages
      rather than lowess.

\item Changes to binGLMDispersion() to allow trended dispersion for data
      sets with small numbers of genes, but with extra warnings.

\item dispDeviance() and dispPearson() now give graceful estimates and
      messages when the dispersion is outside the specified interval.

\item Bug fix to mglmOneWay(), which was confusing parametrizations when
      the design matrix included negative values.

\item mglmOneWay() (and hence glmFit) no longer produces NA coefficients
      when some of the fitted values were exactly zero.

\item Changes to offset behaviour in estimateGLMCommonDisp(),
      estimateGLMTrendedDisp() and estimateGLMTagwiseDisp() to fix bug.
      Changes to several other functions on the way to fixing bugs
      when computing dispersions in data sets with genes that have all
      zero counts.

\item Bug fix to mglmSimple() with matrix offset.
    
\item Bug fix to adjustedProfLik() when there are fitted values exactly
      at zero for one or more groups.
}}

