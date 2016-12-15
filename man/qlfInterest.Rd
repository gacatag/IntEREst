\name{qlfInterest}
\alias{qlfInterest}
\title{
	quasi-likelihood F-test
}
\description{
	Compute quasi-likelihood F-test using \code{edgeR} package. For more information see \code{glmQLFit} and \code{\link[edgeR]{glmQLFTest}} functions in \code{edgeR} package.
}
\usage{
	qlfInterest(x, design=c(), silent=TRUE, disp="common", 
		coef=c(), contrast=NULL, 
		poisson.bound=TRUE, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{design}{
	Design matrix.
}
\item{silent}{
	Whether run silently, i.e. without printing the top differential expression tags. The default is TRUE.
}
\item{disp}{
	The method of estimating the dispersion in the data. Available options are: "common", "trended", "tagwiseInitCommon" and "tagwiseInitTrended". It is also possible to assign a number.
}
\item{coef}{
	Integer or character vector indicating which coefficients of the linear model are to be tested equal to zero. See \code{\link[edgeR]{glmQLFTest}} for more information.
}
\item{contrast}{
	Numeric vector or matrix specifying contrasts of the linear model coefficients to be tested equal to zero. See \code{\link[edgeR]{glmQLFTest}} for more information.
}
\item{poisson.bound}{
	Logical value, if TRUE (i.e. default) the pvalue would be higher than when obtained fom likelihood ratio test while Negative Binomial dispersion is zero. 
}
\item{...}{
	Other parameter settings for the \code{\link[edgeR]{glmQLFTest}} function in the \code{edgeR} package.
}

}
\value{
All values produced by \code{\link[edgeR]{glmQLFTest}} plus the following :
\item{dispersionType}{The name of the type of dispersion used.}
\item{dispersion}{The estimated dispersion values.}
}
\seealso{
	\code{\link{exactTestInterest}}, \code{\link{glmInterest}}, \code{\link{treatInterest}}
}

\author{Ali Oghabian}

\examples{

#Test retention differentiation across the 3 types of sampels
group <- getAnnotation(mdsChr22Obj)[,"type"]
qlfRes<- qlfInterest(x=mdsChr22Obj, 
	design=model.matrix(~group), silent=TRUE, 
	disp="tagwiseInitTrended", coef=2:3, contrast=NULL)

}

