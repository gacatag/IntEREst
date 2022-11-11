\name{subInterestResult}
\alias{subInterestResult}
\title{
Extract subset of object 
}
\description{
Build a new object using subset of data in an \code{SummarizedExperiment} 
object. 
}
\usage{
subInterestResult(x, selectRow, selectCol, 
	sampleAnnoCol, sampleAnnotation=c())
}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{selectRow}{
Numeric or TRUE/FALSE Vector indicating what rows to extract.
}
\item{selectCol}{
A vector with Numeric values, character strings (sample names) or TRUE/FALSE 
Vector indicating what columns to extract.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of object \code{x} to consider for subset data 
extraction.
}
\item{sampleAnnotation}{
Vector including the annotations to consider for subset data extraction. They 
should be present in the \code{sampleAnnoCol} column of the \code{colData} of 
\code{x}.
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
readFreqColIndex<- grep("_readCnt$",colnames(interestDat))
scaledRetentionColIndex<- grep("_fpkm$",colnames(interestDat))
samNames<-paste("sam", 1:4, sep="")
frqTmp<-as.matrix(interestDat[, readFreqColIndex])
sclTmp<-as.matrix(interestDat[, scaledRetentionColIndex])
colnames(frqTmp)<- samNames
colnames(sclTmp)<- samNames
interestResObj<- InterestResult(
	resultFiles=paste("file",1:4, sep="_"),
	rowData= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	counts= frqTmp,
	scaledRetention= sclTmp ,
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy"),
		row.names=samNames
	)
)

#Build new object
newObj<- subInterestResult(interestResObj, selectRow=1:20)

#View newObj
print(newObj)

}

