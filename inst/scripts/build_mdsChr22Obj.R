#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop(
"the path to the bam files and should be supplied as the 1st argument and 
the path to the ouput data should be supplied as the 2nd argument", 
call.=FALSE)
}
MDS_Chr22_BAMFILES<-dir(path = args[1], pattern = ".bam$", full.names = TRUE)
names(MDS_Chr22_BAMFILES)<-
	gsub("\\.bam", "", gsub("/","",gsub(args[1], "",MDS_Chr22_BAMFILES)))

library("IntEREst")

### INTRON RETENTION MEASUREMENT

# Building reference; Do collapse exons
refseqRef<- referencePrepare (sourceBuild="UCSC", ucscGenome="hg19",
    addCollapsedTranscripts=TRUE, ignore.strand=FALSE)

# Annotate U12 type introns
library("BSgenome.Hsapiens.UCSC.hg19")
refseqAnnoMat<- annotateU12(
    pwmU12U2= list(
        pwmU12db[[1]][,11:17],
        pwmU12db[[2]],
        pwmU12db[[3]][,38:40],
        pwmU12db[[4]][,11:17],
        pwmU12db[[5]][,38:40]),
    pwmSsIndex= list(
        indexDonU12=1, 
        indexBpU12=1, 
        indexAccU12=3, 
        indexDonU2=1, 
        indexAccU2=3), 
    referenceChr= refseqRef[,"chr"], 
    referenceBegin= as.numeric(refseqRef[,"begin"]), 
    referenceEnd= as.numeric(refseqRef[,"end"]), 
    referenceIntronExon= refseqRef[,"int_ex"],
    intronExon= "intron",
    matchWindowRelativeUpstreamPos= c(NA,-29,NA,NA,NA),
    matchWindowRelativeDownstreamPos= c(NA,-9,NA,NA,NA), 
    minMatchScore= c( rep("80%", 2), "40%", "80%",  "40%"), 
    refGenome= BSgenome.Hsapiens.UCSC.hg19, 
    setNaAs= "U2", 
    annotateU12Subtype= TRUE)


# Creating temp directory to store the results
outDir<- file.path(tempdir(),"interest_sequential_results")
dir.create(outDir)
outDir<- normalizePath(outDir)

# Summarize read sequences in bam files. Measure Intron retention levels
for(i in 1:length(MDS_BAMFILES)){
    dir.create(paste(outDir, names(MDS_BAMFILES)[i],
        sep="/"), recursive = TRUE)
   interest(
        bamFileYieldSize=1000000,
        bamFile=MDS_BAMFILES[i],
        isPaired=TRUE,
        isPairedDuplicate=FALSE,
        isSingleReadDuplicate=NA,
        reference=refseqRef,
        referenceGeneNames=refseqRef[,"collapsed_transcripts_id"],
        referenceIntronExon=refseqRef[,"int_ex"],
        repeatsTableToFilter= c(),
        junctionReadsOnly=FALSE,
        outFile=paste(outDir, names(MDS_BAMFILES)[i], 
            "interestRes.tsv", sep="/"),
        logFile=paste(outDir, names(MDS_BAMFILES)[i], 
            "log.txt", sep="/"),
        method="IntRet",
        clusterNo=40,
        returnObj=FALSE, 
        scaleLength= TRUE, 
        scaleFragment= TRUE
    )
}

# Reading the intron retention results
# and build SummarizedExperiment object
mdsRefObj<-readInterestResults(
    resultFiles=paste(outDir, names(MDS_BAMFILES), 
            "interestRes.tsv", sep="/"), 
    sampleNames=names(MDS_BAMFILES), 
    sampleAnnotation=data.frame( 
        type=c(rep("ZRSR2mut",8), rep("ZRSR2wt",4), rep("HEALTHY",4)),
        test_ctrl=c(rep("test",8), rep("ctrl",8))), 
    commonColumns=1:9, freqCol=10, scaledRetentionCol=11,
    scaleLength=TRUE, scaleFragment=TRUE, reScale=TRUE, 
    geneIdCol="collapsed_transcripts_id")

# update the object with the intron type (U12- or U2-type) annotations
mdsRefObj<- updateRowDataCol(mdsRefObj,  "intron_type", refseqAnnoMat[,1])

# Choose U12 gebnes on chromosome 22 only
rowInd<- which(
	rowData(mdsRefObj)$collapsed_transcripts_id %in%
		unique(rowData(mdsRefObj)$collapsed_transcripts_id[
			which(rowData(mdsRefObj)$intron_type=="U12" & 
				rowData(mdsRefObj)$chr=="chr22")]))
mdsChr22Obj<-mdsRefObj[rowInd,]

save(mdsChr22Obj,file=paste(args[2],"mdsChr22Obj.rda",sep="/"))




### EXON-EXON JUNCTION MEASUREMENT
library("IntEREst")
library("BSgenome.Hsapiens.UCSC.hg19")
refseqUncollapsed<- referencePrepare (sourceBuild="UCSC", ucscGenome="hg19",
    addCollapsedTranscripts=FALSE, collapseExons=FALSE, ignore.strand=FALSE)

# Union exons of transcripts with overlapping exons
# Keep one copy from each repeating set of exons, only
refExDf<- unionRefTr(referenceChr= refseqUncollapsed[,"chr"], 
    referenceBegin= as.numeric(refseqUncollapsed[,"begin"]), 
    referenceEnd= as.numeric(refseqUncollapsed[,"end"]), 
    referenceTr=as.character(refseqUncollapsed[,"transcript_id"]),
    referenceIntronExon=refseqUncollapsed[,"int_ex"],
    intronExon="exon",
    silent=FALSE)


# Creating temp directory to store the results
outDir<- file.path(tempdir(),"interest_sequential_results")
dir.create(outDir)
outDir<- normalizePath(outDir)

# Summarize read sequences in bam files. Measure exon-exon junction levels
for(i in 1:length(MDS_BAMFILES)){
    dir.create(paste(outDir, names(MDS_BAMFILES)[i],
        sep="/"), recursive = TRUE)
   interest(
        bamFileYieldSize=1000000,
        bamFile=MDS_BAMFILES[i],
        isPaired=TRUE,
        isPairedDuplicate=FALSE,
        isSingleReadDuplicate=NA,
        reference=refExDf,
        referenceGeneNames=refExDf[,"transcripts_id"],
        referenceIntronExon=refExDf[,"int_ex"],
        repeatsTableToFilter= c(),
        junctionReadsOnly=TRUE,
        outFile=paste(outDir, names(MDS_BAMFILES)[i], 
            "exExRes.tsv", sep="/"),
        logFile=paste(outDir, names(MDS_BAMFILES)[i], 
            "log.txt", sep="/"),
        method=c("ExEx"),
        clusterNo=20,
        returnObj=FALSE, 
        scaleLength= FALSE, 
        scaleFragment= TRUE,
        bpparam= SnowParam(workers=25)
    )
}

# Reading the intron retention results
# and build SummarizedExperiment object
mdsExRefObj<-readInterestResults(
    resultFiles=paste(outDir, names(MDS_BAMFILES), 
            "exExRes.tsv", sep="/"), 
    sampleNames=names(MDS_BAMFILES), 
    sampleAnnotation=data.frame( 
        type=c(rep("ZRSR2mut",8), rep("ZRSR2wt",4), rep("HEALTHY",4)),
        test_ctrl=c(rep("test",8), rep("ctrl",8))), 
    commonColumns=1:5, freqCol=6, scaledRetentionCol=7,
    scaleLength=FALSE, scaleFragment=TRUE, reScale=FALSE, 
    geneIdCol="transcripts_id")

# Choose U12 genes on chromosome 22 only
trChooseGr<- GRanges(
	tapply(rowData(mdsChr22Obj)$chr,rowData(mdsChr22Obj)$
		collapsed_transcripts_id,unique),
    IRanges(tapply(rowData(mdsChr22Obj)$begin,rowData(mdsChr22Obj)$
		collapsed_transcripts_id,min),
        tapply(rowData(mdsChr22Obj)$end,rowData(mdsChr22Obj)$
			collapsed_transcripts_id,max)))

mdsExRefObjGr<- GRanges(rowData(mdsExRefObj)$chr,
    IRanges(rowData(mdsExRefObj)$begin,
        rowData(mdsExRefObj)$end))

exOverlap<- findOverlaps(mdsExRefObjGr, trChooseGr, type="any")


mdsChr22ExObj<- mdsExRefObj[unique(queryHits(exOverlap)),]
save(mdsChr22ExObj,file=paste(args[2],"mdsChr22ExObj.rda",sep="/"))

