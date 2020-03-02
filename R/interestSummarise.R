# Function to correct length  based on possible removal of repeat elements
correctLengthRepeat<-function(ref, repeatsTableToFilter){	

	lenRef=ref[,"end"]-ref[,"begin"]+1
	if(length(repeatsTableToFilter)!=0){
		repeatGRanges= GenomicRanges::GRanges( seqnames=
			repeatsTableToFilter[,"chr"],
			IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
			end=repeatsTableToFilter[,"end"], 
			width=repeatsTableToFilter[,"end"]-
				repeatsTableToFilter[,"begin"]+1))
		reduceRepeatGr=GenomicRanges::reduce(repeatGRanges)
		reduceRepeatBegin=as.numeric(GenomicRanges::start(reduceRepeatGr))
		reduceRepeatBeginGr= GenomicRanges::GRanges(
			seqnames=as.character(GenomicRanges::seqnames(reduceRepeatGr)),
				IRanges::IRanges(start=reduceRepeatBegin,
					end=reduceRepeatBegin))

		reduceRepeatEnd=as.numeric(GenomicRanges::end(reduceRepeatGr))
		reduceRepeatEndGr= GenomicRanges::GRanges(
			seqnames=as.character(GenomicRanges::seqnames(reduceRepeatGr)), 
			 	IRanges::IRanges(start=reduceRepeatEnd, end=reduceRepeatEnd))

		refGr=GenomicRanges::GRanges( seqnames=ref[,"chr"],
			IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"],
				width=ref[,"end"]-ref[,"begin"]+1))
		#Repeat elements that are fully located in the reference (Intron/Exon)
		hitsRepeat= GenomicRanges::findOverlaps(reduceRepeatGr, refGr, 
			type= "within")
		#Repeat elements which their 'end' coordinates are only located in the
		#reference
		hitsEndRepeat= GenomicRanges::findOverlaps(reduceRepeatEndGr, refGr, 
			type= "within")
		hitsEndRepeatFilQuery= S4Vectors::queryHits(hitsEndRepeat)[
			is.na(match(S4Vectors::queryHits(hitsEndRepeat), 
				S4Vectors::queryHits(hitsRepeat) ))]
		hitsEndRepeatFilSubject=S4Vectors::subjectHits(hitsEndRepeat)[
			is.na(match(S4Vectors::queryHits(hitsEndRepeat), 
				S4Vectors::queryHits(hitsRepeat) ))]
		#Repeat elements which their 'begin' coordinate are only located in the
		#reference
		hitsBeginRepeat= GenomicRanges::findOverlaps(reduceRepeatBeginGr, 
			refGr, type= "within")
		hitsBeginRepeatFilQuery= S4Vectors::queryHits(hitsBeginRepeat)[
			is.na(match(S4Vectors::queryHits(hitsBeginRepeat), 
				S4Vectors::queryHits(hitsRepeat) ))]
		hitsBeginRepeatFilSubject= S4Vectors::subjectHits(hitsBeginRepeat)[
			is.na(match(S4Vectors::queryHits(hitsBeginRepeat), 
				S4Vectors::queryHits(hitsRepeat) ))]

		lenRef[S4Vectors::subjectHits(hitsRepeat)]= 
			lenRef[S4Vectors::subjectHits(hitsRepeat)]- 
				as.numeric(GenomicRanges::width(reduceRepeatGr)
					[S4Vectors::queryHits(hitsRepeat)])
		if(length(hitsEndRepeatFilQuery)>0){
			lenRef[hitsEndRepeatFilSubject]=
				lenRef[hitsEndRepeatFilSubject]-
					(reduceRepeatEnd[hitsEndRepeatFilQuery]-
						ref[hitsEndRepeatFilSubject,"begin"]+1)
		}
		if(length(hitsBeginRepeatFilQuery)>0){
			lenRef[hitsBeginRepeatFilSubject]=lenRef[hitsBeginRepeatFilSubject]-
				(ref[hitsBeginRepeatFilSubject,"end"]-
					reduceRepeatBegin[hitsBeginRepeatFilQuery]+1)
		}

	}
	return(lenRef)
}

# Function for summarizing the results
interestSummarise <-function(
	reference,
	referenceIntronExon,
	inAnRes,
	method=c("IntRet","ExEx", "IntSpan"),
	referenceGeneNames,
	outFile,
	logFile="",
	appendLogFile=TRUE,
	repeatsTableToFilter=c(),
	scaleLength= c(TRUE,TRUE), 
	scaleFragment= c(TRUE,TRUE))
{
	if(logFile!="")
		cat( "IntERESt:interestSummarise: Begins ...\n", file=logFile, 
			append=appendLogFile)
	cat( "IntERESt:interestSummarise: Begins ...\n")

	if(is(reference,"GRanges")){
		if(length(names(reference))>0){
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=as.numeric(GenomicRanges::end(reference)), 
				strand=as.character(GenomicRanges::strand(reference)), 
				names=as.character(names(reference)))
		} else{
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=as.numeric(GenomicRanges::end(reference)), 
				strand=as.character(GenomicRanges::strand(reference)))
		}
		reference=tmpReference
	}

	res=reference
	writeResults=FALSE
	intExInd=which(method=="IntRet")
	exExInd=which(method=="ExEx")
	intSpanInd=which(method=="IntSpan")

	inAnMat<- matrix(inAnRes, 
		ncol=length(which(c(length(intExInd)>0, length(exExInd)>0, 
			length(intSpanInd)>0))), 
			byrow=FALSE)

	if(length(intExInd)>0){
		ref=reference[which(referenceIntronExon=="intron"),]

		# Calculate the corrected length
		
		lenRef=correctLengthRepeat(ref, repeatsTableToFilter)
		writeResults=TRUE
		res<- cbind(res, inAnMat[ , intExInd])
		colnames(res)[ncol(res)]="IntRet_frequency"
		msg<-
"IntERESt:interestSummarise: Normalizing intron retention read levels.\n"
		if(logFile!="")
			cat( msg, file=logFile, append=TRUE)
		cat(msg)
		frq<- inAnMat[ which(referenceIntronExon=="intron"), intExInd]
		referenceGeneNamesTmp<- as.character(referenceGeneNames[
			referenceIntronExon=="intron"])
		geneCnt<- tapply(frq,referenceGeneNamesTmp,sum)
		FPKM=frq
		if(scaleFragment[intExInd])
			FPKM=((10^6)*frq)/(as.numeric(geneCnt[referenceGeneNamesTmp])+1)
		if(scaleLength[intExInd])
			FPKM=((10^3)*FPKM)/lenRef

		resTmp=rep(0,nrow(reference))
		resTmp[referenceIntronExon=="intron"]=as.numeric(FPKM)
		res=cbind(res, resTmp)
		colnames(res)[ncol(res)]="IntRet_genewide_scaled"
	} 
	if(length(exExInd)>0) {
		ref=reference[referenceIntronExon=="exon",]
		# Calculate the corrected length
		lenRef=correctLengthRepeat(ref, repeatsTableToFilter)

		writeResults=TRUE
		res<- cbind(res, inAnMat[ , exExInd])
		colnames(res)[ncol(res)]="ExEx_frequency"
		msg<-
"IntERESt:interestSummarise: Normalizing exon-exon junction read levels.\n"
		if(logFile!="")
			cat( msg, file=logFile, append=TRUE)
		cat( msg)
		frq<- inAnMat[ which(referenceIntronExon=="exon"), exExInd]
		referenceGeneNamesTmp=as.character(referenceGeneNames[
			which(referenceIntronExon=="exon")])
		geneCnt=tapply(frq,referenceGeneNamesTmp,sum)
		FPKM=frq
		if(scaleFragment[exExInd])
			FPKM=((10^6)*FPKM)/(as.numeric(geneCnt[referenceGeneNamesTmp])+1)
		if(scaleLength[exExInd])
			FPKM=((10^3)*FPKM)/lenRef

		resTmp=rep(0,nrow(reference))
		resTmp[referenceIntronExon=="exon"]=as.numeric(FPKM)
		res=cbind(res, resTmp)
		colnames(res)[ncol(res)]="ExEx_genewide_scaled"
	}
	if(length(intSpanInd)>0){
		ref=reference[which(referenceIntronExon=="intron"),]

		# Calculate the corrected length
		
		lenRef=correctLengthRepeat(ref, repeatsTableToFilter)
		writeResults=TRUE
		res<- cbind(res, inAnMat[ , intSpanInd])
		colnames(res)[ncol(res)]="IntSpan_frequency"
		msg<-
"IntERESt:interestSummarise: Normalizing intron retention read levels.\n"
		if(logFile!="")
			cat( msg, file=logFile, append=TRUE)
		cat(msg)
		frq<- inAnMat[ which(referenceIntronExon=="intron"), intSpanInd]
		referenceGeneNamesTmp<- as.character(referenceGeneNames[
			referenceIntronExon=="intron"])
		geneCnt<- tapply(frq,referenceGeneNamesTmp,sum)
		FPKM=frq
		if(scaleFragment[intSpanInd])
			FPKM=((10^6)*frq)/(as.numeric(geneCnt[referenceGeneNamesTmp])+1)
		if(scaleLength[intSpanInd])
			FPKM=((10^3)*FPKM)/lenRef

		resTmp=rep(0,nrow(reference))
		resTmp[referenceIntronExon=="intron"]=as.numeric(FPKM)
		res=cbind(res, resTmp)
		colnames(res)[ncol(res)]="IntSpan_genewide_scaled"
	} 



	if(writeResults){
		write.table(res, outFile, col.names=TRUE, row.names=FALSE, sep='\t', 
			quote=FALSE)
	} else {
		stop(
'wrong method setting. The correct values are "IntRet", "ExEx", and "IntSpan".'
			)
	}
}
