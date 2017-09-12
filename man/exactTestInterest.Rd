\name{exactTestInterest}
\alias{exactTestInterest}
\title{
	Exact test
}
\description{
Compute genewise exact test between two groups of read counts, using the 
\code{edgeR} package. 
}
\usage{
exactTestInterest(x, sampleAnnoCol=c(), sampleAnnotation=c(), 
	geneIdCol, silent=TRUE, group=c(), rejection.region="doubletail", 
	big.count=900, prior.count=0.125, disp="common", ...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of \code{x} to consider for the analysis.
}
\item{sampleAnnotation}{
A vector of size 2 which cotains values from \code{colData} of 
\code{SummarizedExperiment} object; e.g. if \code{getAnnotation(x)[,
sampleAnnoCol]= c("test", "test", "ctrl","ctrl", ...)} , and the goal is to 
compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be 
\code{c("test", "ctrl")} or \code{c("ctrl", "test")}.
}
\item{geneIdCol}{
Column name (or number of column) in \code{rowData} of \code{x}, i.e. 
\code{SummarizedExperiment} object, that represents the gene ID of the 
introns and exons in \code{x}.
}
\item{silent}{
Whether run the function silently, i.e. without printing the top differential 
expression tags.
}
\item{group}{
Vector to manually define the sample groups (or annotations). It is ignored if 
\code{sampleAnnopCol} is defined.
}
\item{rejection.region}{
The \code{rejection.region} parameter in \code{\link[edgeR]{exactTest}} from 
\code{edgeR} package.
}
\item{big.count}{
The \code{big.count} parameter in \code{\link[edgeR]{exactTest}} from 
\code{edgeR} package.
}
\item{prior.count}{
The \code{prior.count} parameter in \code{\link[edgeR]{exactTest}} from 
\code{edgeR} package.
}
\item{disp}{
The type of estimating the dispersion in the data. Available options are: 
"tagwise", "trended", "common" and "genewise". It is also possible to assign a 
number for manually setting the \code{disp}.
}
\item{...}{
Other parameter settings for the \code{\link[edgeR]{estimateDisp}} function 
(e.g. the \code{design} parameter) in the \code{edgeR} package.
}

}
\value{
\item{table}{ 
Data frame containing columns for the log2 fold-change (logFC), the average of 
log2 counts-per-million (logCPM), and the two-sided p-value (PValue).  
}
\item{comparison}{The name of the two compared groups.}
\item{dispersionType}{The name of the type of dispersion used.}
\item{dispersion}{The estimated dispersion values.}
}
\seealso{
\code{\link{lfc}}, \code{\link{glmInterest}}, \code{\link{qlfInterest}}, 
\code{\link{treatInterest}}, \code{\link{DEXSeqIntEREst}}
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

res<- exactTestInterest(InterestResultObj, sampleAnnoCol="gender", 
	sampleAnnotation=c("F","M"), geneIdCol= "gene_id", 
	silent=TRUE, disp="common")

}

