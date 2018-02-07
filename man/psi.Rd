\name{psi}
\alias{psi}
\title{
	Psi values estimation
}
\description{
Calculating Psi values (percentage spliced in) base on two count matrices from 
a single or two separate objects.  
}
\usage{
psi (x, y, intCol, exCol, pseudoCnt=0)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{y}{ 
	Optional; i.e. an object of type \code{SummarizedExperiment}.
}

\item{intCol}{
Column numbers or column names in counts matrix of \code{x} which include the 
number of reads mapped to the introns.}
\item{exCol}{
Column numbers or column names in counts matrix of \code{x} (or if defined 
\code{y}) which include the number of reads spanning the introns (or mapping
exons flanking the introns).
}
\item{pseudoCnt}{Pseudo counts to sum to the denominator of the devision to 
avoid devision to zero.
}

}
\value{
data.frame with column size equal to the size of \code{intCol} parameter, and
row size equal to the number of rows in \code{x}. It contains the psi values 
(i,e.values between 0 and 1 showing the fraction of spliced in transcripts).
}

\seealso{
	\code{\link{interestResultIntEx}}
}

\author{Ali Oghabian}

\examples{

mdsChr22IntObj<- mdsChr22Obj[which(rowData(mdsChr22Obj)$int_ex=="intron"), ]

#Build object including intron-retention and exon-junction results
mdsChr22RefIntExObj<- interestResultIntEx(intObj=mdsChr22Obj, 
	exObj=mdsChr22ExObj, mean.na.rm=TRUE, postExName="ex_junc", 
	intExCol="int_ex" )
# Calculate Psi
psiRes<- psi(mdsChr22RefIntExObj,
	intCol=which(colData(mdsChr22RefIntExObj)$intronExon=="intron"),
	exCol=which(colData(mdsChr22RefIntExObj)$intronExon=="exon"))
# show Psi results
head(psiRes)
}

