\name{InterestResult}
\alias{InterestResult}
\title{
	Building SummarizedExperiment object from results in IntEREst.
}
\description{
Calls the constructors and creates a \code{SummarizedExperiment} object. 
For more information on the resulted object and the class see 
\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
}
\usage{

InterestResult(resultFiles=c(), counts, scaledRetention, 
	scaleLength, scaleFragment, sampleAnnotation, rowData)

}
\arguments{
\item{resultFiles}{
Vector of link to the result files of \code{interest}.
}
\item{counts}{
Numeric Matrix that includes the read counts.
}
\item{scaledRetention}{
Matrix that includes the scaled retention values.
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
\item{rowData}{ 
Data frame with Intron/Exon annotations and read count  and scaled retention 
values for each sample.
}
}
\value{
Returns an object of class \code{SummarizedExperiment}. 
}
\seealso{
	\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
	\code{\link{attributes}}
	\code{\link{addAnnotation}}
	\code{\link{counts-method}}
	\code{\link{plot-method}}
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
		gender=c("M","M","F","F"), row.names=paste("sam", 1:4, sep="")
	)
)

# View object
InterestResultObj

}

