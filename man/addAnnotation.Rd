\name{addAnnotation}
\alias{addAnnotation}
\title{
	Adding sample annotations to a InterestResult object
}
\description{
	Adds a new sample annotation to the InterestResult object.
}
\usage{
addAnnotation(x, sampleAnnotationType, sampleAnnotation)

}
\arguments{
\item{x}{ 
	Object of type \code{InterestResult}.
}
\item{sampleAnnotationType}{
	The name of the new column to be added to the \code{@sampleAnnotation} table of \code{InterestResult} object.
}
\item{sampleAnnotation}{
	Vector with the same length as the row-size of the \code{@sampleAnnotation} attribute of the \code{InterestResult} object, which includes the sample annotations.
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

