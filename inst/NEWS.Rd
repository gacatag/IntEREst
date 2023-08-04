\name{IntEREstnews}
\title{IntEREst News}
\encoding{UTF-8}


\section{CHANGES IN VERSION 1.26.0}{
\subsection{SIGNIFICANT USER-VISIBLE CHANGES}{

		\itemize{

			\item \code{IntSpan} method is changed to count reads that span the whole 
			intron precisely. Previous these were counted separayted for 5' end and 
			3' end of the introns.
			
			\item \code{repeatsTableToFilter} is ignored for IntSpana nd ExSkip as 
			the skipped regions of the reads (denoted by N or sometimes D in the 
			CIGAR of bam/sam file ) are taken into count for these benchmarks (not 
			the mapped regions which are  dentoed with Ms in the cigar).

		}
	}
}

\section{CHANGES IN VERSION 1.24.0}{

        \subsection{NEW FEATURES}{

                \itemize{

                        \item \code{strandSpecific} is a new parameter added to \code{interest()} and \code{interest.sequential()}. It indicates that
IntEREst is now strand-aware. All analysis can be run whilst taking into account (or ignoring) the strand specificity of the sequencing reads. 

                        \item \code{referencePrepare} can create references (with collapsed exons) that does not ignore the strand information.

                }
        }

	\subsection{SIGNIFICANT USER-VISIBLE CHANGES}{

		\itemize{

			\item \code{GenomicFiles::reduceByYield} function is now used instead of 
			\code{BiocParallel::bpiterat} for improved efficiency and readability. 
			This relates to how analysis related to different pieces of the alignment 
			data is distributed over the parallel cores and the results are sumarized.

		}
	}
	\subsection{BUG FIXES}{

      \itemize{

              \item In \code{interest()} and and \code{interest.sequential()} 
              The oredr of the results of the different methods in the result 
              matrix is corrected. Correct running of analyses using multiple 
              methods with a single command is now possible.

      }
 }
}

\section{CHANGES IN VERSION 1.22.0}{

	\subsection{NEW FEATURES}{

		\itemize{
			
			\item \code{method} parameter in \code{interest()} and \code{interest.sequential()} functions now support 
\code{"ExSkip"}. For each exon, the fucntions can now count the number of reads that skip over the exons.

		}
	}

}

\section{CHANGES IN VERSION 1.20.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{limitRanges} parameter in \code{interest()} and \code{interest.sequential()} functions allow for targetted analysis. It only loads reads that map to the defined coordinates.
			
		}
	}

}

\section{CHANGES IN VERSION 1.12.0}{

        \subsection{BUG FIXES}{

                \itemize{

                        \item Replacing \code{class( x ) == "foo"} with \code{is( x , "foo")} .

                }
        }
}

\section{CHANGES IN VERSION 1.10.0}{

	\subsection{NEW FEATURES}{

		\itemize{

			\item \code{applyOverlap()} Apply a function on counts of a 
			\code{SummarizedExperiment} basd on the overlaps of its genome 
			coordinates to with another \code{SummarizedExperiment}.
		}
	}

	\subsection{BUG FIXES}{

		\itemize{

			\item \code{referencePrepare()} stopped from modifying chromosome 
			names.I doesn't add "chr" to the beginning of chromosome names.

		}
	}
}

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




