\name{IntEREstnews}
\title{IntEREst News}
\encoding{UTF-8}


\section{CHANGES IN VERSION 1.2.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{DEXSeqIntEREst()} runs DEXSeq differential exon usage
			or intron retention test on the SummarizedExperiment objects 
			(resulted from \code{interest()} or \code{interest.sequential()}
			analysis).

			\item \code{unionRefTr()} performs union on the genes in reference
			data frame with overlapping introns/exons. The resulted data frame 
			features from each repeating exon or intron a single copy only.

			\item \code{deseqInterest()} runs differential intron retention 
			test adapted from the \code{DESeq2} package.

			\item \code{interestResultIntEx()} builds a 
			\code{SummarizedExperiment} object from an intron retention and an
			exon-exon junction object (resulted from \code{interest()}, 
			\code{interest.sequential()} and/or \code{InterestResult()} 
			functions). The average of the junction levels of the flanking 
			exons are added to the \code{SummerizedExperiment} with the intron 
			retention values.

		}
	}

	\subsection{SIGNIFICANT USER-VISIBLE CHANGES}{

		\itemize{

			\item \code{u12Index()} supports \code{intronExon} parameter that 
			provides the possibility to extract either the rows that represent
			the U12-type introns or the exons flanking to the U12-type introns
			from a \code{SummerizedExperiment} object.

		}
	}

	\subsection{BUG FIXES}{

		\itemize{

			\item Correcting interestAnalyse.R and interestAnalyse.sequential.R
			so that \code{interest()} and \code{interest.sequential()} support 
			bam files with single (unpaired) reads.

			\item Correct bpparam setting in  \code{interest()}.

		}
	}
}
