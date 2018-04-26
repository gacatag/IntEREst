\name{IntEREstnews}
\title{IntEREst News}
\encoding{UTF-8}

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

		\item \code{interest.sequential()} and \code{interest()} corrections
		to their object output option.

		\item \code{annotateU12()} modified to work correctly with the new 
		changes in \code{Biostrings} package.

		}
	}
}




