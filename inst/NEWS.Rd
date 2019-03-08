\name{IntEREstnews}
\title{IntEREst News}
\encoding{UTF-8}


\section{CHANGES IN VERSION 1.8.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{interest()} and \code{interest.sequential()} functions 
			support the possibility to define other parameter settings of 
			\code{BamFile()} function of \code{Rsamtools} package e.g. 
			"qnamePrefixEnd" and "qnameSuffixStart".

			\item \code{annotateU12()} can return the PWM match scores (if 
			requested). The hybrid AT-AG and GT-AC U12 types can be set to be
			ignored. The introns/exons coordinates can be filtered according 
			to their map to another set of coordinates.

		}
	}

	\subsection{BUG FIXES}{

		\itemize{

			\item \code{referencePrepare()} just gives warning and continue
			without annotating gene names (when \code{annotateGeneIds=TRUE}
			parameter is set) if annotating with gene names is not
			possible.

			\item \code{interest()} used to return Error if single end mapped 
			reads did not existed at all (despite being requetsed). This has 
			been corrected.

		}
	}
}

\section{CHANGES IN VERSION 1.6.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{buildSsTypePwms()} function supports the possibility 
			to select begin and end point of splice sites sequences from which
			PWMs are built. It also supports the possibility to paste splice 
			sites to build PWM.

		}
	}

	\subsection{BUG FIXES}{

		\itemize{

			\item \code{interest.sequential()} and \code{interest()}: 
			corrections	to their object output option.
	
			\item \code{annotateU12()} modified to work correctly with the new 
			changes in \code{Biostrings} package.

			\item \code{buildSsTypePwms()} corrected and modified to better 
			suit data for all species from SpliceRack and U12DB.

		}
	}
}

\section{CHANGES IN VERSION 1.4.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{interest()} and \code{interest.sequential()} functions 
			now support \code{"IntSpan"} method, allows counting intron-
			spanning reads.

			\item \code{psi()} function is added. It calculates Psi values.

			\item \code{annotateU12()} fucntion supports DNAStringSet objects 
			as its \code{refGenome} input.

			\item Improved vignette document.


		}
	}

	\subsection{SIGNIFICANT USER-VISIBLE CHANGES}{

		\itemize{

			\item \code{eBayesInterest()} removed.

		}
	}

	\subsection{BUG FIXES}{

		\itemize{

			\item \code{interest.sequential()} and \code{interest()}: 
			corrections to their object output option.

			\item \code{annotateU12()} modified to work correctly with the new
			changes in \code{Biostrings} package.

		}
	}
}




