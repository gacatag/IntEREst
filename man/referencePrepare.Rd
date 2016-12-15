\name{referencePrepare}
\alias{referencePrepare}
\title{Creates reference file}
\description{Creates reference file for IntEREst functions, e.g. \code{interest()}. The function uses functions of \code{biomaRt} library.}
\usage{
referencePrepare( outFileTranscriptsAnnotation="",
	annotateGeneIds=TRUE, 
	u12IntronsChr=c(), u12IntronsBeg=c(), u12IntronsEnd=c(),
	u12IntronsRef,	collapseExons=TRUE, sourceBuild="UCSC", 
	ucscGenome="hg19", ucscTableName="knownGene",
	ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/",
	biomart="ENSEMBL_MART_ENSEMBL",
	biomartDataset="hsapiens_gene_ensembl",
	biomartTranscriptIds=NULL, biomartExtraFilters=NULL, 
	biomartIdPrefix="ensembl_",	biomartHost="www.ensembl.org",
	biomartPort=80,	circSeqs="", miRBaseBuild=NA, taxonomyId=NA,
	filePath="", fileFormat=c("auto", "gff3", "gtf"), fileDatSrc=NA,
	fileOrganism=NA, fileChrInf=NULL, 
	fileDbXrefTag=c(), addCollapsedTranscripts=TRUE, 
	ignore.strand=FALSE )
}
\arguments{
\item{outFileTranscriptsAnnotation}{If defined outputs transcripts annotations.}
\item{annotateGeneIds}{Wether annotate and add the gene ids information.}
\item{collapseExons}{Whether collapse (i.e. reduce) the exonic regions. TRUE by default.}
\item{sourceBuild}{The source to use to build the reference data, \code{"UCSC"}, \code{"biomaRt"}, and \code{"file"} (for GFF3 or GTF files) are supported.}
\item{ucscGenome}{The genome to use. \code{"hg19"} is the default. See \code{genome} parameter of \code{makeTxDbFromUCSC} function of \code{edgeR} library for more information.}
\item{ucscTableName}{The UCSC table name to use. See \code{tablename} parameter of \code{makeTxDbFromUCSC} function of \code{edgeR} library for more information.}
\item{ucscUrl}{The UCSC URL address.  See \code{url} parameter of \code{makeTxDbFromUCSC} function of \code{edgeR} library for more information.}
\item{u12IntronsChr}{A vector of character strings that includes chromsomal locations of the U12 type introns. If defined together with \code{u12IntronsBeg} and \code{u12IntronsBeg}, they would be used to annotate the U12-type introns.}
\item{u12IntronsBeg}{A vector of numbers that defines the begin (or start) coordinates of the u12-type introns.}
\item{u12IntronsEnd}{A vector of numbers that defines the end coordinates of the u12-type introns.}
\item{u12IntronsRef}{A GRanges object that includes the coordinates of the U12 type introns. If defined, it would be used to annotate the U12-type introns.}
\item{biomart}{BioMart database name. See \code{biomart} parameter of \code{makeTxDbFromBiomart} function of \code{edgeR} library for more information.}
\item{biomartDataset}{BioMart dataset name; default is "hsapiens_gene_ensembl". See \code{dataset} parameter of \code{makeTxDbFromBiomart} function of \code{edgeR} library for more information.}
\item{biomartTranscriptIds}{optional parameter to only retrieve transcript annotation results for a defined set of transcript ids. See \code{transcript_ids} parameter of \code{makeTxDbFromBiomart} function of \code{edgeR} library for more information.}
\item{biomartExtraFilters}{A list of names; i.e. additional filters to use in the BioMart query. See \code{filters} parameter of \code{makeTxDbFromBiomart} function of \code{edgeR} library for more information.}
\item{biomartIdPrefix}{A list of names; i.e. additional filters to use in the BioMart query. See \code{id_prefix} parameter of \code{makeTxDbFromBiomart} function of \code{edgeR} library for more information.}
\item{biomartHost}{ Host to connect to; the default is "www.ensembl.org". For older versions of the GRCH you can provide the archive websites, e.g. for GRCH37 you can use "grch37.ensembl.org".}
\item{biomartPort}{The port to use in the HTTP communication with the host. Default is 80.}
\item{circSeqs}{A character vector that includes chromosomes that should be marked as circular. See \code{circ_seqs} parameter of \code{makeTxDbFromBiomart} and \code{makeTxDbFromUCSC} functions of \code{edgeR} library for more information.}
\item{miRBaseBuild}{Set appropriate build Information from mirbase.db to use for microRNAs (default=NA). See \code{miRBaseBuild} parameter of \code{makeTxDbFromBiomart} and \code{makeTxDbFromUCSC} functions of \code{edgeR} library for more information.}
\item{taxonomyId}{This parameter can be used to provide taxonomy Ids. It is set to Na by default. You can check the taxonomy Ids with the \code{available.species()} function in \code{GenomeInfoDb} package. For more information see \code{taxonomyId} parameter of \code{makeTxDbFromBiomart} and \code{makeTxDbFromUCSC} functions of \code{edgeR} library.}
\item{filePath}{Character string i.e. the path to file. Used if \code{sourceBuild} is \code{"file"}. }
\item{fileFormat}{The format of the input file. \code{"auto"}, \code{"gff3"} and \code{"gtf"} is supported. }
\item{fileDatSrc}{Character string describing the source of the data file. Used if \code{sourceBuild} is \code{"file"}.}
\item{fileOrganism}{The genus and species name of the organism. Used if \code{sourceBuild} is \code{"file"}.}
\item{fileChrInf}{Dataframe that includes information about the chromosome. The first column represents the chromosome name and the second column is the length of the chromosome. Used if \code{sourceBuild} is \code{"file"}.}
\item{fileDbXrefTag}{A vector of chracater strings which if defined it would be used as feature names. Used if \code{sourceBuild} is \code{"file"}.}
\item{addCollapsedTranscripts}{Whether add a column that includes the collapsed transcripts information. Used if \code{collapseExons} is \code{TRUE}.}
\item{ignore.strand}{Whether consider the strands in the reference. If set \code{TURE} the strands would be ingnored.}
}
\value{
	Data frame that includes the coordinates and annotations of the introns and exons of the transcripts, i.e. the reference.
}

\author{Ali Oghabian}
\examples{

	# Build test gff3 data
	tmpGen<- u12[u12[,"ens_trans_id"]=="ENST00000413811",]
	tmpEx<-tmpGen[tmpGen[,"int_ex"]=="exon",]
	exonDat<- cbind(tmpEx[,3], ".", 
		tmpEx[,c(7,4,5)], ".", tmpEx[,6], ".",paste("ID=exon", 
		tmpEx[,11], "; Parent=ENST00000413811", sep="") )
	trDat<- c(tmpEx[1,3], ".", "mRNA", as.numeric(min(tmpEx[,4])), 
		as.numeric(max(tmpEx[,5])), ".", tmpEx[1,6], ".", 
		"ID=ENST00000413811")

	outDir<- file.path(tempdir(),"tmpFolder")
	dir.create(outDir)
	outDir<- normalizePath(outDir)

	gff3File=paste(outDir, "gffFile.gff", sep="/")

	cat("##gff-version 3\n",file=gff3File, append=FALSE)
	cat(paste(paste(trDat, collapse="\t"),"\n", sep=""),
		file=gff3File, append=TRUE)

	write.table(exonDat, gff3File,
		row.names=FALSE, col.names=FALSE,
		sep='\t', quote=FALSE, append=TRUE)	

	# Selecting U12 introns info from 'u12' data
	u12Int<-u12[u12$int_ex=="intron"&u12$int_type=="U12",]

	# Test the function
	refseqRef<- referencePrepare (sourceBuild="file", 
		filePath=gff3File, u12IntronsChr=u12Int[,"chr"], 
		u12IntronsBeg=u12Int[,"begin"], 
		u12IntronsEnd=u12Int[,"end"], collapseExons=TRUE, 
		fileFormat="gff3", annotateGeneIds=FALSE)
}
