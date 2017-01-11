\name{attributes}
\alias{attributes}
\alias{getAnnotation}
\alias{scaledRetention}
\alias{nread}
\alias{interestDf}

\title{
 Extracting values of useful attributes of \code{InterestResult} objects
}
\description{
Several functions are provided that can extract various attribnutes from an 
object of class \code{InterestResult}. It is possible to extract sample 
annotations using \code{getAnnotation} function. One can also extract the 
retention levels scaled to the length of the introns/exons, fragments mapped to
the genes or both (dependant on user setting) using \code{scaledRetention} 
function. The \code{nread} function returns the read counts and 
\code{interestDf} extracts intron/exon annotations.
}
\usage{
getAnnotation(x)

scaledRetention(x)
nread(x)
interestDf(x)
}
\arguments{
\item{x}{ 
	Object of type \code{InterestResult}.
}
}
\value{
Various data types (data.frame/vector) dependent on the function used. See the 
"Description" for more information.
}
\seealso{
	\code{\link{addAnnotation}}
}

\author{Ali Oghabian}

\examples{

# Retrieve the sample annotations from mdsChr22Obj
getAnnotation(mdsChr22Obj)
# Retrieving the read counts from mdsChr22Obj
head(nread(mdsChr22Obj))
# Retrieving the scaled retention levels from mdsChr22Obj
head(scaledRetention(mdsChr22Obj))
# Retrieving the annotations and read counts from mdsChr22Obj
head(interestDf(mdsChr22Obj))

}

