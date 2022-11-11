\name{pwmU12db}
\alias{pwmU12db}
\docType{data}
\title{
	PWM of U12 and U2-type introns splice sites
}
\description{
PWM of U12 and U2-type introns splice sites and it is based on the U12DB 
database. 
}
\usage{data("pwmU12db")}
\format{
A list that contains Position Weight Matrices (PWM) of donor site, branch 
point and acceptor site of U12-type introns and the PWMs of donor site and 
acceptor site of U2-type introns. It is based on the U12DB database.
  \describe{
    \item{\code{pwmDonU12}}{
		A position weigh matrix for the donor site of the U12-type introns, 
		with 4 rows and 46 columns. The rows of the matrix represent "A", "C", 
		"G", and "T" nucleotides and the columns represent the postions in the 
		genome. Each position in the matrix include a weight (i.e. number 
		between 0 and 1) which indicates how common the corresponding base 
		(represented by the row of the matrix) is observed in the corerespoding
		position (represented by the colum of the matrix).}
    \item{\code{pwmBpU12}}{
		A position weigh matrix for the branch point of the U12-type introns, 
		with 4 rows and 9 columns.}
    \item{\code{pwmAccU12}}{
		A position weigh matrix for the acceptor site of the U12-type introns, 
		with 4 rows and 46 columns.}
    \item{\code{pwmDonU2}}{
		A position weigh matrix for the donor site of the U2-type introns, with
		4 rows and 25 columns.}
    \item{\code{pwmAccU2}}{
		A position weigh matrix for the acceptor site of the U12-type introns, 
		with 4 rows and 46 columns.}
  }
}
\value{
List of 5 numeric matrices representing the PWMs of donor site of U12-type 
introns, branch point site of U12-type introns, acceptor site of U12-type 
introns, donor site of U2-type introns, and acceptor site of U2-type introns.
}
\source{
Alioto, T.S. U12DB: a database of orthologous U12-type spliceosomal introns. 
Nucleic Acids Research 2006, doi: 10.1093/nar/gkl796
}
\keyword{datasets}
