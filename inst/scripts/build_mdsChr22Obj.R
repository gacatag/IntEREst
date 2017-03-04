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

# Creating temp directory to store the results
outDir<- file.path(tempdir(),"interest_sequential_results")
dir.create(outDir)
outDir<- normalizePath(outDir)

# Analyze the bam files
lapply(1:length(MDS_Chr22_BAMFILES), function(i){

    dir.create(paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
        sep="/"))
    interest.sequential(
        bamFileYieldSize=10000,
        bamFile=MDS_Chr22_BAMFILES[i],
        isPaired=TRUE,
        isPairedDuplicate=FALSE,
        isSingleReadDuplicate=NA,
        reference=u12,
        referenceGeneNames=u12[,"ens_gene_id"],
        referenceIntronExon=u12[,"int_ex"],
        repeatsTableToFilter=c(),
        outFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
            "interestRes.tsv", sep="/"),
        logFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
            "log.txt", sep="/"),
        method=c("IntRet","ExEx"),
        returnObj=FALSE,
        scaleLength= c(TRUE,FALSE), 
        scaleFragment= c(TRUE,TRUE)
    )
    }
)

#Reading the intron retention results
mdsChr22Obj<-readInterestResults(
    resultFiles=paste(outDir, names(MDS_Chr22_BAMFILES), 
            "interestRes.tsv", sep="/"), 
    sampleNames=names(MDS_Chr22_BAMFILES), 
    sampleAnnotation=data.frame( 
        type=c(rep("ZRSR2_Mut",8), rep("ZRSR2_WT",4), rep("Healthy",4)),
        test_ctrl=c(rep("test",8), rep("ctrl",8))), 
    commonColumns=1:16, freqCol=17, scaledRetentionCol=18,
    scaleLength=TRUE, scaleFragment=TRUE, reScale=FALSE)

# Choose chromosome 22 only
index<- which(rowData(mdsChr22Obj)[,"chr"]=="chr22")
mdsChr22ObjFil<- subInterestResult(mdsChr22Obj, selectRow=index)

mdsChr22Obj<- mdsChr22ObjFil

save(mdsChr22Obj,file=paste(args[2],"mdsChr22Obj.rda",sep="/"))

