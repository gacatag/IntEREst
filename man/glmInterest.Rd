\name{glmInterest}
\alias{glmInterest}
\title{
	generalized linear model likelihood ratio tests
}
\description{
	Compute generalized linear model likelihood ratio tests using \code{edgeR} package. For more information see \code{\link[edgeR]{glmfit}} and \code{glmLRT()} functions in \code{edgeR} package.
}
\usage{
	glmInterest(x, design=c(), silent=TRUE, disp="common", 
		coef=c(), contrast=NULL, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{design}{
	Design matrix.
}
\item{silent}{
	Whether run silently, i.e. without printing the top differential expression tags. Default is TRUE.
}
\item{disp}{
	The method of estimating the dispersion in the data. Available options are: "common", "trended", "tagwiseInitCommon" and "tagwiseInitTrended". It is also possible to assign a number.
}
\item{coef}{
	Integer or character vector indicating which coefficients of the linear model are to be tested equal to zero. See \code{glmLRT()} in \code{edgeR} for more information.
}
\item{contrast}{
	Numeric vector or matrix specifying contrasts of the linear model coefficients to be tested equal to zero. See \code{glmLRT()} in \code{edgeR} for more information.
}
\item{...}{
	Other parameter settings for the \code{glmLRT()} function in the \code{edgeR} package.
}

}
\value{
All values produced by \code{glmLRT} in \code{edgeR} package plus the following :
\item{dispersionType}{The name of the type of dispersion used.}
\item{dispersion}{The estimated dispersion values.}
}
\seealso{
	\code{\link{exactTestInterest}}, \code{\link{qlfInterest}}, \code{\link{treatInterest}}
}

\author{Ali Oghabian}

\examples{

#Test retention differentiation across the 3 types of sampels
group <- getAnnotation(mdsChr22Obj)[,"type"]
glmRes<- glmInterest(x=mdsChr22Obj, 
	design=model.matrix(~group), silent=TRUE, 
	disp="tagwiseInitTrended", coef=2:3, contrast=NULL)

}

