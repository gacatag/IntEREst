\name{mergeInterestResult}
\alias{mergeInterestResult}
\title{
	merge two SummarizedExperiment objects into one 
}
\description{
	Build a new object bu merging data of two \code{SummarizedExperiment} objects.
}
\usage{
mergeInterestResult(x, y)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{y}{ 
	Object of type \code{SummarizedExperiment}.
}
}
\value{
An object of calss \code{SummarizedExperiment}.
}

\seealso{
	\code{\link{interest}},
	\code{\link{InterestResult}}.
}

\author{Ali Oghabian}

\examples{


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
		sam1_fpkm=fpkm1,
		sam2_fpkm=fpkm2,
		sam3_fpkm=fpkm3,
		sam4_fpkm=fpkm4
)
readFreqColIndex<- grep("_readCnt$",colnames(interestDat))
scaledRetentionColIndex<- grep("_fpkm$",colnames(interestDat))

scalRetTmp<- as.matrix(interestDat[ ,scaledRetentionColIndex])
colnames(scalRetTmp)<-gsub("_fpkm$","", colnames(scalRetTmp))

frqTmp<- as.matrix(interestDat[ ,readFreqColIndex])
colnames(frqTmp)<-gsub("_readCnt$","", colnames(frqTmp))

#Object including data for Males
interestResObjM<-InterestResult(
	resultFiles=paste("file",1:2, sep="_"),
	rowData= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	counts= frqTmp[,1:2],
	scaledRetention= scalRetTmp[,1:2],
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:2, sep=""),
		gender=c("M","M"),
		health=c("healthy","unhealthy"),
		row.names=paste("sam", 1:2, sep="")
	)
)

#Object including data for Females
interestResObjF<-InterestResult(
	resultFiles=paste("file",3:4, sep="_"),
	rowData= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	counts= frqTmp[,3:4],
	scaledRetention= scalRetTmp[,3:4],
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",3:4, sep=""),
		gender=c("F","F"),
		health=c("healthy","unhealthy"),
		row.names=paste("sam", 3:4, sep="")
	)
)

#Build new object
newObj<- mergeInterestResult(interestResObjM, interestResObjF)

#View newObj
print(newObj)

}

