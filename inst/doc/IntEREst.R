## ----pipeline, out.width = 600, ppi=350, fig.retina = NULL, fig.align="center", echo=FALSE, eval=TRUE, fig.cap="**Figure 1:** Diagram of IntEREst running pipeline"----
knitr::include_graphics("../inst/fig/IntEREst.png")

## ----parallel_pipeline, out.width = 500, echo=TRUE, eval=FALSE-----------
#  library("IntEREst")
#  library("MDS.Chr22.U12Genes")
#  
#  # Creating temp directory to store the results
#  outDir<- file.path(tempdir(),"interest_parallel_results")
#  dir.create(outDir)
#  outDir<- normalizePath(outDir)
#  
#  # Analyze the bam files
#  lapply(1:length(MDS_Chr22_BAMFILES), function(i){
#  	dir.create(paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  		sep="/"))
#  	interest(
#  		bamFileYieldSize=1000000,
#  		tmpDir=paste(outDir, names(MDS_Chr22_BAMFILES)[i], "tmp",
#  			sep="/"),
#  		bamFile=MDS_Chr22_BAMFILES[i],
#  		filterPairedDuplicate=TRUE,
#  		filterSingleReadDuplicate=FALSE,
#  		reference=u12,
#  		referenceGeneNames=u12[,"ens_gene_id"],
#  		referenceIntronExon=u12[,"int_ex"],
#  		repeatsTableToFilter=c(),
#  		outFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  			"interestRes.tsv", sep="/"),
#  		logFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  			"log.txt", sep="/"),
#  		delTmpFolder=TRUE,
#  		method=c("IntRet","ExEx"),
#  		clusterNo=3,
#  		returnObj=FALSE,
#  		scaleLength= c(TRUE,FALSE),
#  		scaleFragment= c(TRUE,TRUE)
#  	)
#  	}
#  )
#  

## ----sequential_pipeline, out.width = 500, echo=TRUE, eval=FALSE---------
#  library("IntEREst")
#  library("MDS.Chr22.U12Genes")
#  
#  # Creating temp directory to store the results
#  outDir<- file.path(tempdir(),"interest_sequential_results")
#  dir.create(outDir)
#  outDir<- normalizePath(outDir)
#  
#  # Analyze the bam files
#  lapply(1:length(MDS_Chr22_BAMFILES), function(i){
#  
#  	dir.create(paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  		sep="/"))
#  	interest.sequential(
#  		bamFileYieldSize=1000000,
#  		tmpDir=paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  			"tmp", sep="/"),
#  		bamFile=MDS_Chr22_BAMFILES[i],
#  		filterPairedDuplicate=TRUE,
#  		filterSingleReadDuplicate=FALSE,
#  		reference=u12,
#  		referenceGeneNames=u12[,"ens_gene_id"],
#  		referenceIntronExon=u12[,"int_ex"],
#  		repeatsTableToFilter=c(),
#  		outFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  			"interestRes.tsv", sep="/"),
#  		logFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i],
#  			"log.txt", sep="/"),
#  		delTmpFolder=TRUE,
#  		method=c("IntRet","ExEx"),
#  		returnObj=FALSE,
#  		scaleLength= c(TRUE,FALSE),
#  		scaleFragment= c(TRUE,TRUE)
#  	)
#  	}
#  )

## ----build_object, out.width = 500, echo=TRUE, eval=FALSE----------------
#  #Reading the intron retention results
#  mdsChr22Obj<-readInterestResults(
#  	resultFiles=paste(outDir, names(MDS_Chr22_BAMFILES),
#  			"interestRes.tsv", sep="/"),
#  	sampleNames=names(MDS_Chr22_BAMFILES),
#  	sampleAnnotation=data.frame(
#  		type=c(rep("ZRSR2_Mut",8), rep("ZRSR2_WT",4), rep("Healthy",4)),
#  		test_ctrl=c(rep("test",8), rep("ctrl",8))),
#  	commonColumns=1:16, freqCol=17, scaledRetentionCol=18,
#  	scaleLength=TRUE, scaleFragment=TRUE, reScale=FALSE)
#  
#  # Choose chromosome 22 only
#  index<- which(interestDf(mdsChr22Obj)[,"chr"]=="chr22")
#  mdsChr22ObjFil<- subInterestResult(mdsChr22Obj,
#      interestDfRow=index, interestDfSample=c(), sampleAnnotation=2,
#      sampleAnnoCol=c()
#  )
#  
#  mdsChr22Obj<- mdsChr22ObjFil
#  

## ----view_object, out.width = 500, echo=TRUE, eval=TRUE------------------
# Load library quietly
suppressMessages(library("IntEREst"))
#View object
print(mdsChr22Obj)

## ----plot_intron_object, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center", fig.cap="**Figure 2:** Plotting the distribution of the retention levels ($log_e$ scaled retention) of introns of genes located on chromosome 22. The values have been averaged across the sample types ZRSR2 mutated, ZRSR2 wild type, and healthy."----

# Retention of all introns
plot(mdsChr22Obj, logScaleBase=exp(1), pch=20, loessLwd=1.2, 
	summary="mean", col="black", sampleAnnoCol="type", 
	lowerPlot=TRUE, upperPlot=TRUE)

## ----plot_u12intron_object, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center", fig.cap="**Figure 3:** Plotting the distribution of the retention levels ($log_e$ scaled retention) of introns of genes located on chromosome 22. The values have been averaged across the sample types ZRSR2 mutated, ZRSR2 wild type, and healthy."----

#Retention of U12 introns
plot(mdsChr22Obj, logScaleBase=exp(1), pch=20, plotLoess=FALSE, 
	summary="mean", col="black", sampleAnnoCol="type", 
	subsetRows=u12Index(mdsChr22Obj))

## ----exact_test, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center"----
# Check the sample annotation table
getAnnotation(mdsChr22Obj)

# Run exact test
test<- exactTestInterest(mdsChr22Obj, 
	sampleAnnoCol="test_ctrl", sampleAnnotation=c("ctrl","test"), 
	geneIdCol= "ens_gene_id", silent=TRUE, disp="common")

# Number of stabilized introns (in Chr 22)
sInt<- length(which(test$table[,"PValue"]<0.05 
	& test$table[,"logFC"]>0 & 
	interestDf(mdsChr22Obj)[,"int_ex"]=="intron"))
print(sInt)
# Number of stabilized (significantly retained) U12 type introns
numStU12Int<- length(which(test$table[,"PValue"]<0.05 & 
	test$table[,"logFC"]>0 & 
	interestDf(mdsChr22Obj)[,"int_type"]=="U12" & 
	!is.na(interestDf(mdsChr22Obj)[,"int_type"])))
# Number of U12 introns
numU12Int<- 
	length(which(interestDf(mdsChr22Obj)[,"int_type"]=="U12" & 
	!is.na(interestDf(mdsChr22Obj)[,"int_type"]))) 
# Fraction(%) of stabilized (significantly retained) U12 introns
perStU12Int<- numStU12Int/numU12Int*100
print(perStU12Int)
# Number of stabilized U2 type introns
numStU2Int<- length(which(test$table[,"PValue"]<0.05 & 
	test$table[,"logFC"]>0 & 
	interestDf(mdsChr22Obj)[,"int_type"]=="U2" & 
	!is.na(interestDf(mdsChr22Obj)[,"int_type"])))
# Number of U2 introns
numU2Int<- 
	length(which(interestDf(mdsChr22Obj)[,"int_type"]=="U2" & 
	!is.na(interestDf(mdsChr22Obj)[,"int_type"])))
# Fraction(%) of stabilized U2 introns
perStU2Int<- numStU2Int/numU2Int*100
print(perStU2Int)

## ----glr_test, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center"----

# Extract type of samples
group <- getAnnotation(mdsChr22Obj)[,"type"]
group

# Test retention levels' differentiation across 3 types samples
qlfRes<- qlfInterest(x=mdsChr22Obj, 
	design=model.matrix(~group), silent=TRUE, 
	disp="tagwiseInitTrended", coef=2:3, contrast=NULL, 
	poisson.bound=TRUE)

# Extract index of the introns with significant retention changes
ind= which(qlfRes$table$PValue< 0.05)
# Extract introns with significant retention level changes
variedIntrons= interestDf(mdsChr22Obj)[ind,]

#Show first 5 rows and columns of the result table
print(variedIntrons[1:5,1:5])


## ----boxplot_object, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center", fig.cap="**Figure 4:** Boxplot of the retention levels of U12 introns vs U2 introns, summed over samples with similar annotations *i.e.* ZRSR2 mutation, ZRSR2 wild type, or healthy."----
# boxplot U12 and U2-type introns 
par(mar=c(7,4,2,1))
u12Boxplot(mdsChr22Obj, sampleAnnoCol="type", 
	intExCol="int_ex",	intTypeCol="int_type", intronExon="intron", 
	col=rep(c("orange", "yellow"),3) ,	lasNames=3, 
	outline=FALSE, ylab="FPKM", cex.axis=0.8)

## ----u12BoxplotNb_object, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center", fig.cap="**Figure 5:** Boxplot of retention levels of U12 introns vs their up- and down-stream U2 introns across all samples."----
# boxplot U12-type intron and its up/downstream U2-type introns 
par(mar=c(2,4,1,1))
u12BoxplotNb(mdsChr22Obj, sampleAnnoCol="type", lasNames=1,
	intExCol="int_ex", intTypeCol="int_type", intronExon="intron", 
	boxplotNames=c(), outline=FALSE, plotLegend=TRUE, 
	geneIdCol="ens_gene_id", xLegend="topleft", 
	col=c("pink", "lightblue", "lightyellow"), ylim=c(0,1e+06), 
	ylab="FPKM", cex.axis=0.8)

## ----density_plot, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.cap="**Figure 6:** Density plot of the log fold change of U12-type introns, random U2-type introns and U2 introns (up / down / up and down)stream of U12-type introns. "----
u12DensityPlotIntron(mdsChr22Obj, 
	type= c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), 
	fcType= "edgeR", sampleAnnoCol="test_ctrl", 
	sampleAnnotation=c("ctrl","test"), intExCol="int_ex", 
	intTypeCol="int_type", strandCol= "strand", 
	geneIdCol= "ens_gene_id", naUnstrand=FALSE, col=c(2,3,4,5,6), 
	lty=c(1,2,3,4,5), lwd=1, plotLegend=TRUE, cexLegend=0.7, 
	xLegend="topright", yLegend=NULL, legend=c(), randomSeed=10,
	ylim=c(0,0.6), xlab=expression("log"[2]*" fold change FPKM"))

## ----lfc-sig, echo = TRUE, eval=TRUE, message = FALSE, fig.width=6, fig.height=4, fig.align="center"----
# estimate log fold-change of introns 
# by comparing test samples vs ctrl 
# and don't show warnings !
lfcRes<- lfc(mdsChr22Obj, fcType= "edgeR", 
	sampleAnnoCol="test_ctrl",sampleAnnotation=c("ctrl","test"))

# Build the order vector
ord<- rep(1,length(lfcRes))
ord[u12Index(mdsChr22Obj)]=2

# Median of log fold change of U2 introns (ZRSR2 mut. vs ctrl)
median(lfcRes[ord==1])
# Median of log fold change of U2 introns (ZRSR2 mut. vs ctrl)
median(lfcRes[ord==2])

# Run Jockheere Terpstra's trend test
library(clinfun)
jonckheere.test(lfcRes, ord, alternative = "increasing", 
	nperm=1000)

