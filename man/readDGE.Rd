\name{readDGE}
\alias{readDGE}

\title{Read and Merge a Set of Files Containing DGE Data}

\description{Reads and merges a set of text files containing digital gene expression data.}

\usage{readDGE(files, path=NULL, columns=c(1,2), group=NULL, labels=NULL, ...)}

\arguments{ 
\item{files}{character vector of filenames, or alternatively a data.frame with a column containing the file names of the files containing the libraries of counts and, optionally, columns containing the \code{group} to which each library belongs, descriptions of the other samples and other information.}
\item{path}{character string giving the directory containing the files.
The default is the current working directory.}
\item{columns}{numeric vector stating which two columns contain the tag names and counts, respectively}
\item{group}{vector, or preferably a factor, indicating the experimental group to which each library belongs. If \code{group} is not \code{NULL}, then this argument overrides any group information included in the \code{files} argument.}
\item{labels}{character vector giving short names to associate with the libraries. Defaults to the file names.}
\item{...}{other are passed to \code{read.delim}}
}

\details{
Each file is assumed to contained digital gene expression data for one sample (or library), with transcript identifiers in the first column and counts in the second column. Transcript identifiers are assumed to be unique and not repeated in any one file. 
By default, the files are assumed to be tab-delimited and to contain column headings.
The function forms the union of all transcripts and creates one big table with zeros where necessary.
}

\value{DGEList object}

\author{Mark Robinson and Gordon Smyth}

\examples{
#  Read all .txt files from current working directory

\dontrun{files <- dir(pattern="*\\\\.txt$")
RG <- readDGE(files)}
}

\seealso{
\code{\link{DGEList}} provides more information about the \code{DGEList} class and the function \code{DGEList}, which can also be used to construct a \code{DGEList} object, if \code{readDGE} is not required to read in and construct a table of counts from separate files.
}

\keyword{file}
