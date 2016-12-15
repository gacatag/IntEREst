# NOTE! The working directory is the path where this file is located
library("MDS.Chr22")
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
		bamFileYieldSize=1000000,
		tmpDir=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
			"tmp", sep="/"),
		bamFile=MDS_Chr22_BAMFILES[i],
		filterPairedDuplicate=T,
		filterSingleReadDuplicate=F,
		reference=u12,
		referenceGeneNames=u12[,"gene_ens_id"],
		referenceIntronExon=u12[,"int_ex"],
		outFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
			"interestRes.tsv", sep="/"),
		logFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
			"log.txt", sep="/"),
		delTmpFolder=TRUE,
		method=c("IntRet","ExEx"),
		returnObj=FALSE
	)
	}
)

# building results object

mdsChr22Obj<-readInterestResults(
	resultFiles=paste(outDir, names(MDS_Chr22_BAMFILES), 
			"interestRes.tsv", sep="/"), 
	sampleNames=names(MDS_Chr22_BAMFILES), 
	sampleAnnotation=data.frame( 
		type=c(rep("ZRSR2_Mut",8), rep("ZRSR2_WT",4), rep("Healthy",4)),
		test_ctrl=c(rep("test",8), rep("ctrl",8))), 
	commonColumns=1:16, freqCol=17, fpkmCol=18)

save(mdsChr22Obj,file="../../data/mdsChr22Obj.rda")

