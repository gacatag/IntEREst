\name{buildSsTypePwms}
\alias{buildSsTypePwms}
\title{
Building Position Weight Matrices for Splice Sites of U12 and U2 type introns.
}
\description{ 
Builds position Weigh Matrices for the donor and acceptor sites of the U12 and 
U2 type introns, and the branchpoint of the U12 type introns. if 
\code{pdfFileSeqLogos} is defined a pdf is also produced that contains the 
sequence logos of the results. The result is a list that contains PWMs of the 
splice sites of U12 and U2 dependent introns.}
\usage{
buildSsTypePwms( cexSeqLogo=1, pdfWidth=35, pdfHeight=10, tmpDir="./",
	u12dbSpecies="Homo_sapiens", 
	pwmSource="U12DB", 
	u12DonorBegin, u12BranchpointBegin, u12AcceptorBegin, 
	u2DonorBegin, u2AcceptorBegin, u12DonorEnd, 
	u12BranchpointEnd, u12AcceptorEnd, u2DonorEnd, 
	u2AcceptorEnd, pasteSites=FALSE,
	splicerackSsLinks=list(
		U12_AT_AC_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.25", 
		U12_AT_AC_branchpoint=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.26",
		U12_AT_AC_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.29",
		U12_GT_AG_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.22", 
		U12_GT_AG_branchpoint=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.27",
		U12_GT_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.21",
		U2_GC_AG_donor="http://katahdin.mssm.edu/splice/out/9606_logo_file.24",
		U2_GC_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.30", 
		U2_GT_AG_donor="http://katahdin.mssm.edu/splice/out/9606_logo_file.23",
		U2_GT_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.28"),
	u12dbLink="ftp://genome.imim.es/pub/software/u12/u12db_v1_0.sql.gz",
	u12dbDbName="u12db", u12dbDropDb=TRUE,  pdfFileSeqLogos="", 
	removeTempFiles=TRUE, ...)
}
\arguments{
\item{cexSeqLogo}{
Font size of sequence logo plots; used only if  \code{pdfFileSeqLogos} is 
defined.}
\item{pdfWidth, pdfHeight}{
The width and height of the graphics region of the pdf in inches. The default 
values are 35 and 10. }
\item{tmpDir}{Path to directory used for storing temporary files.}
\item{u12dbSpecies}{ 
What species data to use when  getting the data from the U12DB database 
(\code{pwmSource="U12DB"}).}
\item{pwmSource}{ 
The source used to buildSplice Sites of U12 and U2 type introns the PWM for U12
and U2 dependent introns. Default is \code{U12DB}; but also accepts 
\code{SpliceRack}. }
\item{u12DonorBegin, u12DonorEnd}{Integer values. They correspond to the begin 
and end point of the donor sequences of U12-type introns to consider 
(optional).}
\item{u12BranchpointBegin, u12BranchpointEnd}{Integer values. Begin and end 
points of the branch point sequences of U12-type introns (optional).} 
\item{u12AcceptorBegin, u12AcceptorEnd}{Integer values. Begin and end points of
the acceptor sequences of U12-type introns (optional).}
\item{u2DonorBegin, u2DonorEnd}{Integer values. Begin and end points of the 
donor sequences of U2-type introns (optional).} 
\item{u2AcceptorBegin, u2AcceptorEnd}{Integer values. Begin and end points of 
the acceptor sequences of U2-type introns (optional).} 
\item{pasteSites}{Logical. If TRUE the donor, branch point and acceptor seqs 
are pasted before a PWM is built; then the PWMs of each (donor, acceptor and 
bp) are assigned. If FALSE (default) the PWMs for each is built separately.}
\item{splicerackSsLinks}{
A list (or vector) that contains the SpliceRack URL links to the text files 
that contain Position Weigh Matrices of the splice sites of U12 and U2 introns.
This parameter is used only when pwmSource="SpliceRack". You can get the links
to PWM files from this URL (choose logo files with "File" links): 
\url{http://katahdin.mssm.edu/splice/splice_matrix.cgi?database=spliceNew}.
The links should  be defined in the following order: \code{U12_AT_AC_donor}, 
\code{U12_AT_AC_branchpoint}, \code{U12_AT_AC_acceptor},
\code{U12_GT_AG_donor}, \code{U12_GT_AG_branchpoint}, 
\code{U12_GT_AG_acceptor}, \code{U2_GC_AG_donor}, \code{U2_GC_AG_acceptor}, 
\code{U2_GT_AG_donor}, and \code{U2_GT_AG_acceptor}. }
\item{u12dbLink}{
A character string containing the URL for downloading the zipped MySQL dump
file of the U12DB. Used when \code{pwmSource="U12DB"}. }
\item{u12dbDbName}{
Name of the database copy of the U12DB that is build locally. Used when 
\code{pwmSource="U12DB"}.}
\item{u12dbDropDb}{
Drop (or remove) the local copy of the U12DB database at the end of the run.
Used when \code{pwmSource="U12DB"}.}
\item{pdfFileSeqLogos}{
Path to PDF file containing the sequence logos of the results. By default it 
does not produce a file.}
\item{removeTempFiles}{ 
Whether remove temporary files at the end of the run; accepts  \code{TRUE} or  
\code{FALSE} values (default is  \code{TRUE}).}
\item{...}{ 
Authorization arguments needed by the DBMS instance. See the manual for 
\code{dbConnect} of the \code{DBI} package for more info.}
}
\value{
\item{pwmDonorU12}{ 
Matrix (with 4 rows represnting A, C, G, T and \code{n} columns representing
the genomic coordinates) representing the Position Weight Matrix of donor site
of U12-type introns. }
\item{pwmBpU12}{Position Weight Matrix of branchpoint of U12-type introns. }
\item{pwmAccU12}{Position Weight Matrix of acceptor site of U12-type introns. }
\item{pwmDonU2}{Position Weight Matrix of donor site of U2-type introns.}
\item{pwmAccU2}{ Position Weight Matrix of acceptor site of U2-type introns. }
}
\seealso{
  \code{\link{annotateU12}}.
}

\author{Ali Oghabian}

\examples{
# Time demanding function
\dontrun{
#Build temp directory  
tmpDir<- tempdir()


# Creating subdirectory for storing u12db temp files
dir.create(paste(tmpDir, "u12dbTmp", sep="/"))

# Extracting PWMs of Splice Sites of U12 and U2 type introns -
# based on u12db
u12dbPwm<-buildSsTypePwms(
	tmpDir=paste(tmpDir, "u12dbTmp", sep="/"),
	u12dbSpecies="Homo_sapiens",
	resource="U12DB",
	u12dbDbName="u12db",
	u12dbDropDb=TRUE,
	removeTempFiles=TRUE)


# Creating subdirectory for storing SpliceRack temp files
dir.create(paste(tmpDir, "splicerackTmp", sep="/"))

# Extracting PWMs of Splice Sites of U12 and U2 type introns - 
# based on SpliceRack
spliceRackPwm<- buildSsTypePwms(
	tmpDir= paste(tmpDir, "splicerackTmp", sep="/"),
	resource="SpliceRack",
	removeTempFiles=TRUE)
}
}





