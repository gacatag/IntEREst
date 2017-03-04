\docType{methods}
\name{counts-method}
\alias{counts.InterestResults}
\alias{counts,SummarizedExperiment-method}
\alias{counts-method}

\title{
	Counts - method
}
\description{
	Returns the (row) number of reads that are mapped to introns/exons in 
	various samples.
}
\usage{
\S4method{counts}{SummarizedExperiment}(object)
}
\arguments{
\item{object}{ 
	Object of type \code{SummarizedExperiment}.
}
}
\value{
	Returns a numeric matrix.
}
\seealso{
Class:
	\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}


Method:
	 \code{\link{plot-method}}.
}
\author{Ali Oghabian}

\examples{
#Show contents of a InterestResults object included in IntEREst 
head(counts(mdsChr22Obj))

#Make a test InterestResults object
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


InterestResultObj<- InterestResult(
	resultFiles=paste("file",1:4, sep="_"),
	rowData= interestDat[ , -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	counts= frqTmp,
	scaledRetention= scalRetTmp,
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
		, row.names=paste("sam", 1:4, sep="")
	)
)

#Show 
head(counts(InterestResultObj))
}

