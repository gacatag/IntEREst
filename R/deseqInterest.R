deseqInterest<- function(x, design,  
	pAdjustMethod = "BH", sizeFactor=c(), 
	contrast, bpparam, ...){
	if(missing(bpparam)){
		parallel<- FALSE
	} else {
		parallel<-TRUE
	}
	dds<- DESeq2::DESeqDataSetFromMatrix(countData = counts(x), 
		colData = colData(x), design = design)
	if(length(sizeFactor)>0)
		DESeq2::sizeFactors(dds)=sizeFactor
	dds<- DESeq2::DESeq(dds, BPPARAM = bpparam)
	cat("\nResult names that can be used for contrasts are:\n")
	cat(DESeq2::resultsNames(dds))
	cat("\n")
	ddsDiff<- DESeq2::results(dds, 
	contrast=contrast, pAdjustMethod=pAdjustMethod, parallel = parallel, 
	BPPARAM = bpparam, ...)
	return(ddsDiff)
}


























