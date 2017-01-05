\name{mergeInterestResult}
\alias{mergeInterestResult}
\title{
	merge two interestResult objects into one 
}
\description{
	Build a new object bu merging data of two \code{interestResult} objects. 
}
\usage{
mergeInterestResult(x, y)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{y}{ 
	Object of type \code{interestResult}.
}
}
\value{
An object of calss \code{interestResult}.
}

\seealso{
	\code{\link{interest}},
	\code{\link{interestResult}}.
}

\author{Ali Oghabian}

\examples{

geneId<- paste("gene", c(rep(1,7), rep(2,7), rep(3,7), rep(4,7)), 
	sep="_")
readCnt1<- sample(1:100, 28)
readCnt2<- sample(1:100, 28)
readCnt3<- sample(1:100, 28)
readCnt4<- sample(1:100, 28)
fpkm1<- readCnt1/(tapply(readCnt1, geneId, sum))[geneId]
fpkm2<- readCnt2/(tapply(readCnt2, geneId, sum))[geneId]
fpkm3<- readCnt3/(tapply(readCnt3, geneId, sum))[geneId]
fpkm4<- readCnt4/(tapply(readCnt4, geneId, sum))[geneId]

# Creating object using test data
interestDat<-data.frame( 
		int_ex=rep(c(rep(c("exon","intron"),3),"exon"),4),
		int_ex_num= rep(c(1,1,2,2,3,3,4),4),         
		int_type=rep(c(NA,"U2",NA,"U12",NA,"U2",NA),4),
		strand=rep("*",28),
		gene_id= geneId,
		sam1_readCnt=readCnt1,
		sam2_readCnt=readCnt2,
		sam3_readCnt=readCnt3,
		sam4_readCnt=readCnt4,
		sam1_fpkm=fpkm1,
		sam2_fpkm=fpkm2,
		sam3_fpkm=fpkm3,
		sam4_fpkm=fpkm4
)
readFreqColIndex<- grep("_readCnt",colnames(interestDat))
scaledRetentionColIndex<- grep("_fpkm",colnames(interestDat))
#Object including data for Males
interestResObjM<-interestResult(
	resultFiles=paste("file",1:2, sep="_"),
	interestDf= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	readFreq= as.matrix(interestDat[, readFreqColIndex[1:2]]),
	scaledRetention= as.matrix(interestDat[, 
		scaledRetentionColIndex[1:2]]),
	sampleNames=paste("sam",1:2, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:2, sep=""),
		gender=c("M","M"),
		health=c("healthy","unhealthy")
	)
)

#Object including data for Females
interestResObjF<-interestResult(
	resultFiles=paste("file",3:4, sep="_"),
	interestDf= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	readFreq= as.matrix(interestDat[, readFreqColIndex[3:4]]),
	scaledRetention= as.matrix(interestDat[, 
		scaledRetentionColIndex[3:4]]),
	sampleNames=paste("sam",3:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",3:4, sep=""),
		gender=c("F","F"),
		health=c("healthy","unhealthy")
	)
)

#Build new object
newObj<- mergeInterestResult(interestResObjM, interestResObjF)

#View newObj
print(newObj)

}

