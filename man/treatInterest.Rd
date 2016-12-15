\name{treatInterest}
\alias{treatInterest}
\title{
	Differential retention test relative to a threshold
}
\description{
	Compute a genewise statistical test relative to a fold-change threshold using \code{edgeR} package. For more information see \code{\link[edgeR]{glmTreat}} function in \code{edgeR} package.
}
\usage{
	treatInterest(x, design=c(), silent=TRUE, disp="common", 
		coef=c(), contrast=NULL, lfc=0, ...)
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
	Integer or character vector indicating which coefficients of the linear model are to be tested equal to zero. See \code{\link[edgeR]{glmTreat}} for more information.
}
\item{contrast}{
	Numeric vector or matrix specifying contrasts of the linear model coefficients to be tested equal to zero. See \code{\link[edgeR]{glmTreat}} for more information.
}
\item{lfc}{
	Numeric scalar i.e. the log fold change threshold. 
}
\item{...}{
	Other parameter settings for the \code{glmFit} function in the \code{edgeR} package.
}

}
\value{
All values produced by \code{\link[edgeR]{glmTreat}} plus the following :
\item{dispersionType}{The name of the type of dispersion used.}
\item{dispersion}{The estimated dispersion values.}
}
\seealso{
	\code{\link{exactTestInterest}}, \code{\link{qlfInterest}}, \code{\link{glmInterest}}
}

\author{Ali Oghabian}

\examples{

group <- getAnnotation(mdsChr22Obj)[,"type"]

#Test retention differentiation across the 3 types of sampels 
# The log fold change threshold is 0
treatRes<- treatInterest(x=mdsChr22Obj, 
	design=model.matrix(~group), silent=TRUE, 
	disp="tagwiseInitTrended", coef=2:3, contrast=NULL, lfc=0)

}

