applyOverlap<- function(
	query,
	subject,
	type="any",
	replaceValues=FALSE,
	intExCol="int_ex", 
	intronExon="intron",
	sujectGeneNamesCol,
	repeatsTableToFilter=c(),
	scaleFragment=TRUE,
	scaleLength=TRUE,
	unmapValue=0,
	FUN=mean,
	...
)
{
	if(is(query,"SummarizedExperiment")){
		dfQue<- SummarizedExperiment::rowData(query)
		indQue<- which(dfQue[,intExCol]%in%intronExon)
		dfQue<- dfQue[indQue,]
		queGr<- GenomicRanges::GRanges (as.character(dfQue[,1]), 
			IRanges::IRanges(as.numeric(dfQue[,2]), 
			as.numeric(dfQue[,3])))
	}

	dfSub<- SummarizedExperiment::rowData(subject)
	indSub<- which(dfSub[,intExCol]%in%intronExon)
	dfSub<- dfSub[indSub,]
	subGr<- GenomicRanges::GRanges (as.character(dfSub[,1]), 
		IRanges::IRanges(as.numeric(dfSub[,2]), as.numeric(dfSub[,3])))
	
	mapTmp<- GenomicRanges::findOverlaps(queGr, 
			subGr, type= type, select="all", maxgap=0)

	time1<- Sys.time()

	sumAgg<- stats::aggregate(x=counts(query)[indQue[S4Vectors::queryHits(
		mapTmp)],], 
		by=list(indSub[S4Vectors::subjectHits(mapTmp)]), FUN=FUN, ...)
	time2<- Sys.time()
	runTime1<- difftime(time2,time1, units="secs")

	outFreq<- matrix(unmapValue, nrow=nrow(subject), ncol=ncol(subject))
	colnames(outFreq)<- colnames(sumAgg)[-1]
	outFreq[ sumAgg[,1] ,]<- as.matrix(sumAgg[,-1])
	if(replaceValues){
		# Get true length
		trueLen<- correctLen(rowData(subject)[indSub,], 
			rowData(query)[indQue,], repeatsTableToFilter=repeatsTableToFilter)
		lenQue<- rowData(subject)[,"end"]-rowData(subject)[,"begin"]+1
		lenMat<- matrix(rep(lenQue, ncol(subject)), 
			ncol=ncol(subject), byrow=FALSE)
		lenMat[indSub, ]<- matrix(rep(trueLen,ncol(query)), 
			ncol=ncol(query), byrow=FALSE)
		# measure sum of read frequency in transcripts
		out<- subject
		SummarizedExperiment::assays(out)$counts<- outFreq
		sumFrqs<- stats::aggregate(x=counts(out), 
		by=list(rowData(out)[,sujectGeneNamesCol]), sum, na.rm=TRUE)
		tmpMat<- sumFrqs[,-1]
		rownames(tmpMat)<- sumFrqs[,1]
		sumMat<- tmpMat[rowData(out)[,sujectGeneNamesCol],]
		# measure FPKM
		if(scaleFragment)
			FPKM=((10^6)*outFreq)/(sumMat+1)
		if(scaleLength)
			FPKM=((10^3)*FPKM)/lenMat
		rownames(FPKM)<-c()
		SummarizedExperiment::assays(out)$scaledRetention<- 
			as.matrix(FPKM)

	} else {
		out<- outFreq
	}
	return (out)
}


correctLen<-function(ref, coll, repeatsTableToFilter=c()){
	lenRef<- ref[,"end"]-ref[,"begin"]+1
	lenColl<- coll[,"end"]-coll[,"begin"]+1
	lenOut<- lenRef
	if(length(coll)>0){
		collGRanges<- GenomicRanges::GRanges( seqnames=
			coll[,"chr"],
			IRanges::IRanges(start=coll[,"begin"],
			end=coll[,"end"], 
			width=coll[,"end"]-
				coll[,"begin"]+1))

		refGr<- GenomicRanges::GRanges( seqnames=ref[,"chr"],
			IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"],
				width=ref[,"end"]-ref[,"begin"]+1))
		redcollGr<- GenomicRanges::reduce(collGRanges)
		hits <- findOverlaps(redcollGr, refGr, type="within", maxgap=0, 
			ignore.strand=T)
		lenSum<- tapply(lenColl[S4Vectors::queryHits(hits)], 
			S4Vectors::subjectHits(hits), sum)

		if(length(lenSum)>0) {
			lenOut[as.numeric(names(lenSum))]<- as.numeric(lenSum)
		}
	}

	return(lenOut)
}
