\name{InterestResult}
\alias{InterestResult}
\title{
	Results object for interest and interest.sequential functions.
}
\description{
Creates an \code{InterestResult} object that includes information from result 
files created by \code{interest()} and \code{interest.sequential()} functions. 
}
\usage{
InterestResult (resultFiles=c(), readFreq, 
	scaledRetention, sampleNames, scaleLength, 
	scaleFragment, sampleAnnotation, interestDf)
}
\arguments{
\item{resultFiles}{
Vector of link to the result files of \code{interest}.
}
\item{readFreq}{
Numeric Matrix that includes the read counts.
}
\item{scaledRetention}{
Matrix that includes the scaled retention values.
}
\item{sampleNames}{
Vector with same size as \code{resultFiles} which contains the names of the 
samples.
}
\item{scaleLength}{
Logical value, indicating whether the intron/exon retention levels are scaled 
to the length of the introns/exons. 
}
\item{scaleFragment}{
Logical value, indicating whether the intron/exon retention levels are scaled 
to the fragments mapped to the genes. 
}
\item{sampleAnnotation}{
Data frame with the row-size equal to the size of \code{resultFiles} and 
\code{sampleAnnotation}. Each column of the matrix represents annotations for 
the samples. Column name represents annotation name.
}
\item{interestDf}{ 
Data frame with Intron/Exon annotations and read count  and scaled retention 
values for each sample.
}
}
\value{
Returns an object of class \code{InterestResult}. 
}
\seealso{
	\code{\link{attributes}},
	\code{\link{addAnnotation}}.
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
readFreqColIndex=grep("_readCnt$",colnames(interestDat))
scaledRetentionColIndex=grep("_fpkm$",colnames(interestDat))
InterestResultObj<- InterestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDat[,-c(readFreqColIndex, 
		scaledRetentionColIndex)],
	readFreq=as.matrix(interestDat[,readFreqColIndex]),
	scaledRetention=as.matrix(interestDat[,
		scaledRetentionColIndex]),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE, 
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F")
	)
)

# View object
InterestResultObj

}

