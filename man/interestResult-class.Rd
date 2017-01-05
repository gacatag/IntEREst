\docType{class}
\name{interestResult-class}
\alias{interestResult-class}
\title{interestResult Class}
\description{
Class that contains results of the \code{interest()} and \code{interest.sequential()} analysis.
}
\section{Fields}{

\describe{
\item{resultFiles}{Character vector that includes the path to the result files.}

\item{readFreq}{Matrix that includes the number of reads mapped to the studied introns and exons in the samples.}

\item{scaledRetention}{Matrix that includes the scaled retention levels of the introns and exons in the samples.}

\item{sampleNames}{Character vector that includes the name of the samples.}

\item{scaleLength}{
	Logical value, indicating whether the intron/exon retention levels are scaled to the length of the introns/exons. 
}

\item{scaleFragment}{
	Logical value, indicating whether the intron/exon retention levels are scaled to the fragments mapped to the genes. 
}

\item{sampleAnnotation}{Data frame with the row-size equal to the size of \code{@resultFiles} and \code{sampleNames}. Each column of the matrix represents annotations for the samples. Column name represents annotation name.}

\item{interestDf}{The data frame that includes the retention levels of the introns and exons together with their coordinates and annotations.}
}}

\seealso{
Constructors:
\code{\link{readInterestResults}},
\code{\link{interestResult}}.

Methods:
\code{\link{plot-method}}.

Other functions:
\code{\link{plot-method}}.
}

