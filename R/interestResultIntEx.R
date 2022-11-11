interestResultIntEx<- function(intObj, exObj, intExCol=c(), 
	mean.na.rm=TRUE, postExName="ex_junc" ){
	if(length(intExCol)==1){
		intObj<- intObj[which(rowData(intObj)[,intExCol]=="intron"), ]
		exObj<- exObj[which(rowData(exObj)[,intExCol]=="exon"), ]
	}
	if(length(intExCol)==2){
		intObj<- intObj[which(rowData(intObj)[,intExCol[1]]=="intron"), ]
		exObj<- exObj[which(rowData(exObj)[,intExCol[2]]=="exon"), ]
	}
	exEndGr<- GenomicRanges::GRanges(rowData(exObj)$chr,
		IRanges::IRanges(rowData(exObj)$end+1, 
			rowData(exObj)$end+1))
	exBegGr<- GenomicRanges::GRanges(rowData(exObj)$chr,
		IRanges::IRanges(rowData(exObj)$begin-1, 
			rowData(exObj)$begin-1))
	intGr<- GenomicRanges::GRanges(rowData(intObj)$chr,
		IRanges::IRanges(rowData(intObj)$begin, 
			rowData(intObj)$end))

	intExEnd<-GenomicRanges::findOverlaps(intGr, exEndGr, type="start")
	intExBeg<-GenomicRanges::findOverlaps(intGr, exBegGr, type="end")
	exCnt<- counts(exObj)
	meanExVals<- lapply(1:ncol(exCnt), function(j){
		endVals<- tapply(subjectHits(intExEnd), queryHits(intExEnd), 
			function(tmp) mean(exCnt[unique(tmp),j], na.rm=mean.na.rm))
		begVals<- tapply(subjectHits(intExBeg), queryHits(intExBeg), 
			function(tmp) mean(exCnt[unique(tmp),j], na.rm=mean.na.rm))
		return(trunc(apply(cbind(endVals, begVals), 1, mean, 
			na.rm=mean.na.rm)))
	})
	exCnt<- as.data.frame(meanExVals)
	if(all(colnames(counts(intObj)) %in% colnames(counts(exObj))))
		colnames(exCnt)<- paste(colnames(counts(intObj)),postExName,sep="_")

	exSr<- scaledRetention(exObj)
	meanExSrVals<- lapply(1:ncol(exSr), function(j){
		endVals<- tapply(subjectHits(intExEnd), queryHits(intExEnd), 
			function(tmp) mean(exSr[unique(tmp),j], na.rm=mean.na.rm))
		begVals<- tapply(subjectHits(intExBeg), queryHits(intExBeg), 
			function(tmp) mean(exSr[unique(tmp),j], na.rm=mean.na.rm))
		return(apply(cbind(endVals, begVals), 1, mean, na.rm=mean.na.rm))
	})
	exSr<- as.data.frame(meanExSrVals)
	if(all(colnames(counts(intObj)) %in% colnames(counts(exObj))))
		colnames(exSr)<- paste(colnames(counts(intObj)),postExName,sep="_")
	
	annoMat<- rbind(colData(intObj), colData(intObj))
	annoMat<- cbind(as.data.frame(annoMat), 
		c(rep("intron", nrow(colData(intObj))), 
			rep("exon", nrow(colData(intObj)))))
	rownames(annoMat)<- c(colnames(scaledRetention(intObj)), colnames(exCnt))
	colnames(annoMat)[ncol(annoMat)]<- "intronExon"
	resObj<- InterestResult(
		resultFiles=annoMat[,1], 
		counts=cbind(counts(intObj), exCnt), 
		scaledRetention=cbind(scaledRetention(intObj), exSr), 
		scaleLength=intObj@metadata$scaleFragmentscaleLength, 
		scaleFragment=intObj@metadata$scaleFragmentscaleFragment, 
		sampleAnnotation=annoMat[,-1], 
		rowData=rowData(intObj))
	return(resObj)
}
