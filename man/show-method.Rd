\docType{methods}
\name{show-method}
\alias{show,interestResult-method}
\alias{show-method}

\title{
	show - method
}
\description{
	Show \code{interestResult} objects.
}
\usage{
\S4method{show}{interestResult}(object)
}
\arguments{
\item{object}{ 
	Object of type \code{interestResult}.
}
}
\value{
	Returns NULL.
}
\seealso{
Class:
	\code{\link{interestResult-class}}.


Method:
	 \code{\link{plot-method}}.
}
\author{Ali Oghabian}

\examples{
#Show contents of a interestResults object included in IntEREst 
mdsChr22Obj

#Make a test interestResults object
geneId<- paste("gene", c(rep(1,5), rep(2,5), rep(3,5), rep(4,5)), 
	sep="_")
readCnt1<- sample(1:100, 20)
readCnt2<- sample(1:100, 20)
readCnt3<- sample(1:100, 20)
readCnt4<- sample(1:100, 20)
fpkm1<- readCnt1/(tapply(readCnt1, geneId, sum))[geneId]
fpkm2<- readCnt2/(tapply(readCnt2, geneId, sum))[geneId]
fpkm3<- readCnt3/(tapply(readCnt3, geneId, sum))[geneId]
fpkm4<- readCnt4/(tapply(readCnt4, geneId, sum))[geneId]

# Creating object using test data
interestDat<- data.frame( 
		int_ex=rep(c(rep(c("exon","intron"),2),"exon"),4),
		int_ex_num= rep(c(1,1,2,2,3),4),         
		gene_id= geneId,
		sam1_readCnt=readCnt1,
		sam2_readCnt=readCnt2,
		sam3_readCnt=readCnt3,
		sam4_readCnt=readCnt4,
		sam1_fpkm=readCnt1,
		sam2_fpkm=readCnt2,
		sam3_fpkm=readCnt3,
		sam4_fpkm=readCnt4
)

interestResultObj<- interestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDat,
	readFreqColIndex=grep("_readCnt$",colnames(interestDat)),
	scaledRetentionColIndex=grep("_fpkm$",colnames(interestDat)),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE, 
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F")
	)
)

#Show object
interestResultObj
}

