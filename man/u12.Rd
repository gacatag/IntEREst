\name{u12}
\alias{u12}
\docType{data}
\title{
U12 data
}
\description{
Intron/exon annotations of genes featuring U12 introns. It is based on HG19/GRCh37 (converted from hg17/NCBI35). Moreover the u12 genes are based on the U12DB database.
}
\usage{data("u12")}
\format{
  A data frame with 22713 observations on the following 17 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{int_ex_id}}{a character vector}
    \item{\code{chr}}{a character vector}
    \item{\code{begin}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{strand}}{a numeric vector}
    \item{\code{int_ex}}{a character vector}
    \item{\code{trans_type}}{a character vector}
    \item{\code{ens_gene_id}}{a character vector}
    \item{\code{ens_trans_id}}{a character vector}
    \item{\code{int_ex_num}}{a numeric vector}
    \item{\code{gene_name}}{a character vector}
    \item{\code{trans_name}}{a character vector}
    \item{\code{overlap_no}}{a numeric vector}
    \item{\code{int_type}}{a character vector}
    \item{\code{int_subtype}}{a character vector}
  }
}
\value{
	Data frame that includes the coordinates and annotations of the introns and exons of the transcripts, i.e. the reference.
}
\source{
Alioto, T.S. U12DB: a database of orthologous U12-type spliceosomal introns. Nucleic Acids Research 2006, doi: 10.1093/nar/gkl796
}
\keyword{datasets}
