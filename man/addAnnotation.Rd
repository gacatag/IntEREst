\name{addAnnotation}
\alias{addAnnotation}
\title{
	Adding sample annotations to a \code{SummarizedExperiment} object
}
\description{
Adds a new sample annotation to the \code{SummarizedExperiment} object. In 
other words it adds and column with sample annotations to the colData of the 
\code{SummarizedExperiment} object. 
}
\usage{
addAnnotation(x, sampleAnnotationType, sampleAnnotation)

}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{sampleAnnotationType}{
The name of the new column to be added to the \code{colData} table of 
\code{SummarizedExperiment} object.
}
\item{sampleAnnotation}{
Vector with the same length as the row-size of the \code{colData} attribute of 
the \code{SummarizedExperiment} object, which includes the sample annotations.
}

}
\value{
An InterestResult object.
}

\seealso{
	\code{\link{getAnnotation}}
}

\author{Ali Oghabian}

\examples{
# Check the annotation table of mdsChr22Obj data
getAnnotation(mdsChr22Obj)

# Add a new sample annotation
newMdsChr22Obj <- addAnnotation(x=mdsChr22Obj,
	sampleAnnotationType="sample_number",
	sampleAnnotation=1:16
)

# Retrieve annotations of the new object
getAnnotation(newMdsChr22Obj)

}

