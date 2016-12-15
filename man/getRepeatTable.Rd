\name{getRepeatTable}
\alias{getRepeatTable}
\title{Get table of regions with repetetive DNA sequences}
\description{This function returns a data.frame that includes regions with repetetive DNA sequences. These sequences can bias the mapping of the reads to the genome excluding them will remove the bias.}
\usage{
getRepeatTable( dbUser="genome", 
	dbHost="genome-mysql.cse.ucsc.edu",	ucscGenome="hg19",
	ucscTable="rmsk", minLength=0, repFamilyFil="Alu", 
	repFamilyCol="repFamily", repChrCol="genoName", 
	repBegCol="genoStart", repEndCol="genoEnd", 
	repStrandCol="strand", repNameCol="repName", 
	repClassCol="repClass")
}
\arguments{
\item{dbUser}{Database user name; set as "genome" by default.}
\item{dbHost}{Database host address; set as "genome-mysql.cse.ucsc.edu" by default.}
\item{ucscGenome}{The UCSC genome.}
\item{ucscTable}{The UCSC table name. The table with repetetive sequences by default it is set as \code{"rmsk"}.}
\item{minLength}{the minimum length criteria to consider the repetetive sequences. the default setting is \code{0}.}
\item{repFamilyFil}{A vector including the repeats family to consider. By default the \code{"Alu"} elements are considered.}
\item{repFamilyCol}{The name of the column of the input table (\code{ucscTable}) that represents the repeats family.}
\item{repChrCol}{The column (either name or the number of the column) of the input table that represents the Chromosome names.}
\item{repBegCol}{The column of the table that represents the start coordinates.}
\item{repEndCol}{The column of the table that represents the end coordinates.}
\item{repStrandCol}{The column of the table that represents the strand.}
\item{repNameCol}{The column of the table representing the repeats' names.}
\item{repClassCol}{The column of the table representing the repeats' classes.}
}
\value{
	Data frame with columns representing coordinates and annotations of repetitive DNA elements.  
}

\author{Ali Oghabian}

\examples{

# Download table for Alu elemnts in the human genome
	suppressWarnings(repTable<- getRepeatTable(repFamilyFil="Alu", ucscGenome="hg19"))

}

