\name{eBayesInterest}
\alias{eBayesInterest}
\title{
	Limma based statistical test
}
\description{
	Statistical test based on \code{eBayes()} of \code{limma} package to extract introns that are significantly retained in a group of samples (compared to the other samples). The statistical test is performed on the scaled retention levels of the introns/exons which are normalized to the length of the intron/exons and the sum of the retention levels across across the gene.
}
\usage{
	eBayesInterest (x, sampleAnnoCol=c(), sampleAnnotation=c(), 
		group=c(), design=c(), logBase=2, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{sampleAnnoCol}{
	Which colummn of \code{@sampleAnnotation} of object \code{interestResult} to consider for the analysis.
}
\item{sampleAnnotation}{
	A vector of size 2 which cotains values from \code{@sampleAnnotation} of \code{interestResult} object; e.g. if \code{x@sampleAnnotation[,sampleAnnoCol] = c("test", "test", "ctrl","ctrl", ...)} , and the goal is to compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be \code{c("test", "ctrl")} or \code{c("ctrl", "test")}.
}
\item{group}{
	Vector to manually define the sample groups (or annotations). It is ignored if \code{sampleAnnopCol} is defined.
}
\item{design}{
	Design matrix.
}
\item{logBase}{
	A positive or complex number, i.e. the base to which the scaled retention levels are log transformed. If set as NA the values would not be log transformed. The deafult (and recommended) value is 2. 
}
\item{...}{
	Other parameter settings for the \code{\link[limma]{ebayes}} function of the \code{limma} package.
}

}
\value{
	All values produced by \code{\link[limma]{ebayes}} function of the \code{limma} package.
}
\seealso{
	\code{\link{exactTestInterest}}
}

\author{Ali Oghabian}

\examples{

# Intron retention differentiation test (test vs ctrl)
ebayesRes<- eBayesInterest (x=mdsChr22Obj, 
	sampleAnnoCol="test_ctrl", 
	sampleAnnotation=c("ctrl", "test"), logBase=2)

}

