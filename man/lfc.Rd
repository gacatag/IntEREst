\name{lfc}
\alias{lfc}
\title{
	Log fold change
}
\description{
Log fold change estimation and normalized log fold change using \code{edgeR}
 package. 
}
\usage{
lfc(x, fcType="edgeR", sampleAnnoCol=c(), sampleAnnotation=c(), 
	silent=TRUE, group=c(), rejection.region="doubletail", 
	pseudoCnt=1, log2=TRUE, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{fcType}{
Available as "scaledRetention" or "edgeR" (as default) corresponding to either
log fold change of scaled retention values or degeR normalized log fold change
values.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of \code{x} to consider for the analysis.
}
\item{sampleAnnotation}{
A vector of size 2 which cotains values from \code{colData} of 
\code{SummarizedExperiment} object; e.g. if \code{getAnnotation(x)[,
sampleAnnoCol]= c("test", "test", "ctrl","ctrl", ...)} , and the goal is to 
compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be 
\code{c("test", "ctrl")} or \code{c("ctrl", "test")}.
}
\item{silent}{
Whether run \code{\link{exactTestInterest}} silently, without warnings.
}
\item{group}{
Vector to manually define the sample groups (or annotations). It is ignored if
 \code{sampleAnnoCol} is defined.
}
\item{rejection.region}{
The rejection.region parameter in \code{\link[edgeR]{exactTest}}, considered
 only if \code{fcType} is "edgeR".
}
\item{pseudoCnt}{
	Pseudo count for log transformation (default=1).
}
\item{log2}{
Logical value either TRUE (default) or FALSE indicating whether the foldchanges
should be log 2 transformed.
}
\item{...}{
Other parameter settings from the \code{\link{exactTestInterest}} function.
}

}
\value{
Vector including fold change values.
}

\seealso{
	\code{\link{exactTestInterest}}, \code{\link{u12DensityPlotIntron}}
}

\author{Ali Oghabian}

\examples{


lfcFpkm<- lfc(mdsChr22Obj, fcType="scaledRetention",
	sampleAnnoCol="test_ctrl", 
	sampleAnnotation=c("ctrl", "test"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)


lfcEdgeRFpkm<- lfc(mdsChr22Obj, fcType="edgeR",
	sampleAnnoCol="test_ctrl", 
	sampleAnnotation=c("ctrl", "test"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)

}

