interestIntExAnalysePair <- function(
	readTmp,
	reference,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly)
{

	sam1<- as.data.frame(GenomicAlignments::first(readTmp))[,
		c("qname", "rname", "strand", "pos", "qual", "cigar")]
	sam2<- as.data.frame(GenomicAlignments::second(readTmp))[,
		c("qname", "rname", "strand", "pos", "qual", "cigar")]
	colnames(sam1)<- c("qName","rName","strand","pos","qual",
		"cigar")
	colnames(sam2)<- c("qName","rName","strand","pos","qual",
		"cigar")

######
######!!!
	exExRes<- c()
	intRetRes<- c()
	intSpanRes<- c()
	if( length(which(method %in% c("IntRet","ExEx")))>0 ){ 
######!!!
######
		lenRead1= unlist(lapply(sam1$cigar, function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
		}))
		lenRead2= unlist(lapply(sam2$cigar, function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))
		}))
		locRead1=sam1[,"pos"]
		chrRead1=sam1[,"rName"]
		chrRead2=sam2[,"rName"]
		locRead2=sam2[,"pos"]
		frqEventRead1<- lapply(1:length(sam1$cigar), function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam1$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam1$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp]
					,split="")))]
			frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]

	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead1[tmp], locRead1[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return( list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM) )
	
		})
	
		frqEventRead2<- lapply(1:length(sam2$cigar), function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam2$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam2$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam2$cigar[tmp]
				, split="")))]
			frqEventTmp= frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]
	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead2[tmp], locRead2[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return(list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead2[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM))
		})

		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]]))))
		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			as.character(unlist(tmp[[2]])))))
		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]]))))
		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]]))))
		noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[1]]))))
		chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			as.character(unlist(tmp[[2]])))))
		begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[4]]))))
		endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[5]]))))

		noVec=c(noRead1,noRead2)
		chrVec=c(as.character(chrRead1),
			as.character(chrRead2))
		begVec=c(begRead1,begRead2)
		endVec=c(endRead1,endRead2)

		methodNo=which(method=="ExEx")
		if(length(methodNo)>0){
			ref=reference[which(referenceIntronExon=="exon"),]
			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
				IRanges::IRanges(start=begVec, end=endVec))
			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
					width=ref[,"end"]-ref[,"begin"]+1))
			hits1 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "start")
			hits2 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "end")
			hits3 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "equal")
			hitsSubject<- c()
			hitsQuery<- c()
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"], 
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
	
				filtInd1= which(is.na(match(S4Vectors::queryHits(hits1), 
					S4Vectors::queryHits(hitsFilter))))
				filtInd2= which(is.na(match(S4Vectors::queryHits(hits2),  
					S4Vectors::queryHits(hitsFilter))))
				filtInd3= which(is.na(match(S4Vectors::queryHits(hits3),  
					S4Vectors::queryHits(hitsFilter))))
				hits1= hits1[filtInd1,]
				hits2= hits2[filtInd2,]
				hits3= hits3[filtInd3,]
			}
			# Include the reads that completely map to exons if user asked
			if(!junctionReadsOnly){
				hitsExtra= GenomicRanges::findOverlaps(readGRanges, 
					refGRanges, type= "within")
				hitsSubject=S4Vectors::subjectHits(hitsExtra)
				hitsQuery=S4Vectors::queryHits(hitsExtra)
			}
		
			hitsSubject=c(hitsSubject, S4Vectors::subjectHits(hits1),
				S4Vectors::subjectHits(hits2),S4Vectors::subjectHits(hits3))
			hitsQuery=c(hitsQuery, S4Vectors::queryHits(hits1), 
				S4Vectors::queryHits(hits2), S4Vectors::queryHits(hits3))
			refRes= rep(0,nrow(ref))
			if(length(hitsQuery)>0){
				hitsApply= tapply(noVec[hitsQuery], hitsSubject, function(tmp)
					return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			exExRes<-rep(0, nrow(reference))
			exExRes[which(referenceIntronExon=="exon")]<- refRes
	
		}	
		methodNo=which(method=="IntRet")
		if(length(methodNo)>0){

			ref=reference[which(referenceIntronExon=="intron"),]
			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
				IRanges::IRanges(start=begVec, end=endVec))
			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
					width=ref[,"end"]-ref[,"begin"]+1))
			hits <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "any")
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
				filtInd= which(is.na(match(S4Vectors::queryHits(hits), 
					S4Vectors::queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}
			# Filter reads mapped completely to Either introns or exons
			if(junctionReadsOnly){
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					refGRanges, type= "within")
				filtInd=which(is.na(match(S4Vectors::queryHits(hits), 
					S4Vectors::queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}
			refRes= rep(0,nrow(ref))
			if(length(S4Vectors::queryHits(hits))>0){
				hitsApply=tapply(noVec[S4Vectors::queryHits(hits)], 
					S4Vectors::subjectHits(hits), 
					function(tmp) return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			intRetRes<- rep(0,nrow(reference))
			intRetRes[which(referenceIntronExon=="intron")]<- refRes
	
		}
########
########!!!
	}
	methodNo=which(method=="IntSpan")
	if(length(methodNo)>0){

		sam1IndSp<- grep("N",sam1$cigar)
		sam2IndSp<- grep("N",sam2$cigar)

		lenRead1= unlist(lapply(sam1$cigar[sam1IndSp], function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
		}))
		lenRead2= unlist(lapply(sam2$cigar[sam2IndSp], function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))
		}))
		locRead1=sam1[,"pos"]
		chrRead1=sam1[,"rName"]
		chrRead2=sam2[,"rName"]
		locRead2=sam2[,"pos"]
		frqEventRead1<- lapply(sam1IndSp, function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam1$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam1$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp]
					,split="")))]
			frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]
	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead1[tmp], locRead1[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return( list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM) )
	
		})

		frqEventRead2<- lapply(sam2IndSp, function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam2$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam2$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam2$cigar[tmp]
				, split="")))]
			frqEventTmp= frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]
	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead2[tmp], locRead2[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return(list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead2[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM))
		})

		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]]))))
		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			as.character(unlist(tmp[[2]])))))
		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]]))))
		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]]))))
		noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[1]]))))
		chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			as.character(unlist(tmp[[2]])))))
		begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[4]]))))
		endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[5]]))))
		noVec=c(noRead1,noRead2)
		chrVec=c(as.character(chrRead1),
			as.character(chrRead2))
		begVec=c(begRead1,begRead2)
		endVec=c(endRead1,endRead2)



		ref=reference[which(referenceIntronExon=="intron"),]
		readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
			IRanges::IRanges(start=begVec, end=endVec))

		refIntBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
			IRanges::IRanges(start=ref[,"begin"]-1, end=ref[,"begin"]-1))
		refIntEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
			IRanges::IRanges(start=ref[,"end"]+1, end=ref[,"end"]+1))

		hitsBegInt <- GenomicRanges::findOverlaps(readGRanges, 
			refIntBegGRanges, type= "end")
		hitsEndInt <- GenomicRanges::findOverlaps(readGRanges, 
			refIntEndGRanges, type= "start")
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
				filtBegInd= which(is.na(match(S4Vectors::queryHits(hitsBegInt),
					S4Vectors::queryHits(hitsFilter))))
				hitsBegInt=hitsBegInt[filtBegInd,]
				filtEndInd= which(is.na(match(S4Vectors::queryHits(hitsEndInt),
					S4Vectors::queryHits(hitsFilter))))
				hitsEndInt=hitsEndInt[filtEndInd,]
			}
			# Filter reads mapped completely to Either introns or exons

			refRes= rep(0,nrow(ref))
			if(length(c(S4Vectors::queryHits(hitsBegInt), 
				S4Vectors::queryHits(hitsEndInt)))>0){
				hitsApply=tapply(
					noVec[c(S4Vectors::queryHits(hitsBegInt), 
						S4Vectors::queryHits(hitsEndInt))], 
					c(S4Vectors::subjectHits(hitsBegInt), 
						S4Vectors::subjectHits(hitsEndInt)),
					function(tmp) return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			intSpanRes<- rep(0,nrow(reference))
			intSpanRes[which(referenceIntronExon=="intron")]<- refRes
	}

########!!!
########
	finalRes<-c()
	if(length(intRetRes)>0)
		finalRes<-c(finalRes, intRetRes)
	if(length(exExRes)>0){
		finalRes<-c(finalRes, exExRes)
	}
########
########!!!
	if(length(intSpanRes)>0){
		finalRes<-c(finalRes, intSpanRes)
	}
########!!!
########
	return(finalRes)
#End of defining function for paired-mapped reads that runs on parallel cores 
} 
	

# For single reads
interestIntExAnalyseSingle <- function(
	readTmp,
	reference,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly)
{

# Filtering readTmp based on number of mapping targets for each read
	lenTmp<- unlist(sapply(readTmp, length))
	readTmpFil<- readTmp[which(lenTmp<= maxNoMappedReads)]

#Extract information of the mapped pieces of the reads
	sam1<- as.data.frame(readTmpFil)[, c("qname", "rname", "strand", 
		"pos", "qual", "cigar")]
	colnames(sam1)<- c("qName","rName","strand","pos","qual",
		"cigar")
	sam2= NA

######
######!!!
	exExRes<- c()
	intRetRes<- c()
	intSpanRes<- c()
	if( length(which(method %in% c("IntRet","ExEx")))>0 ){ 
######!!!
######
		lenRead1= unlist(lapply(sam1$cigar, function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
		}))

		lenRead2= 0
		locRead1= sam1[,"pos"]
		chrRead1= sam1[,"rName"]
		chrRead2= 0
		locRead2= 0


		frqEventRead1 =lapply(1:length(sam1$cigar), function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam1$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam1$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp]
					, split="")))]
			frqEventTmp= frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]

	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead1[tmp], locRead1[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return( list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM) )

		})

		frqEventRead2=NA
		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]]))))
		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[2]]))))
		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]]))))
		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]]))))
		noRead2= c()
		chrRead2= c()
		begRead2= c()
		endRead2= c()

# analyzing reads
		noVec=c(noRead1,noRead2)
		chrVec=c(as.character(chrRead1),
			as.character(chrRead2))
		begVec=c(begRead1,begRead2)
		endVec=c(endRead1,endRead2)
		exExRes<- c()
		intRetRes<- c()
		methodNo=which(method=="ExEx")
		if(length(methodNo)>0){
			ref=reference[which(referenceIntronExon=="exon"),]
			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
				IRanges::IRanges(start=begVec, end=endVec))
			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
					width=ref[,"end"]-ref[,"begin"]+1))
			hits1 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "start")
			hits2 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "end")
			hits3 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "equal")
			hitsSubject<- c()
			hitsQuery<- c()
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"], 
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
	
				filtInd1= which(is.na(match(S4Vectors::queryHits(hits1), 
					S4Vectors::queryHits(hitsFilter))))
				filtInd2= which(is.na(match(S4Vectors::queryHits(hits2),  
					S4Vectors::queryHits(hitsFilter))))
				filtInd3= which(is.na(match(S4Vectors::queryHits(hits3),  
					S4Vectors::queryHits(hitsFilter))))
				hits1= hits1[filtInd1,]
				hits2= hits2[filtInd2,]
				hits3= hits3[filtInd3,]
			}
			# Include the reads that completely map to exons if user asked
			if(!junctionReadsOnly){
				hitsExtra= GenomicRanges::findOverlaps(readGRanges, 
					refGRanges, type= "within")
				hitsSubject=S4Vectors::subjectHits(hitsExtra)
				hitsQuery=S4Vectors::queryHits(hitsExtra)
			}
		
			hitsSubject=c(hitsSubject, S4Vectors::subjectHits(hits1),
				S4Vectors::subjectHits(hits2),S4Vectors::subjectHits(hits3))
			hitsQuery=c(hitsQuery, S4Vectors::queryHits(hits1), 
				S4Vectors::queryHits(hits2), S4Vectors::queryHits(hits3))
			refRes= rep(0,nrow(ref))
			if(length(hitsQuery)>0){
				hitsApply= tapply(noVec[hitsQuery], hitsSubject, function(tmp)
					return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			exExRes<-rep(0, nrow(reference))
			exExRes[which(referenceIntronExon=="exon")]<- refRes
	
		}
		methodNo=which(method=="IntRet")
		if(length(methodNo)>0){
	
			ref=reference[which(referenceIntronExon=="intron"),]
			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
				IRanges::IRanges(start=begVec, end=endVec))
			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
				IRanges(start=ref[,"begin"], end=ref[,"end"], 
					width=ref[,"end"]-ref[,"begin"]+1))
			hits <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
				type= "any")
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
				filtInd= which(is.na(match(S4Vectors::queryHits(hits), 
					S4Vectors::queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}
			# Filter reads mapped completely to Either introns or exons
			if(junctionReadsOnly){
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					refGRanges, type= "within")
				filtInd=which(is.na(match(S4Vectors::queryHits(hits), 
					S4Vectors::queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}
			refRes= rep(0,nrow(ref))
			if(length(S4Vectors::queryHits(hits))>0){
				hitsApply=tapply(noVec[S4Vectors::queryHits(hits)], 
					S4Vectors::subjectHits(hits), 
					function(tmp) return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			intRetRes<- rep(0,nrow(reference))
			intRetRes[which(referenceIntronExon=="intron")]<- refRes
	
		}
########
########!!!
	}
	methodNo=which(method=="IntSpan")
	if(length(methodNo)>0){

		sam1IndSp<- grep("N",sam1$cigar)
		sam2IndSp<- c()

		lenRead1= unlist(lapply(sam1$cigar[sam1IndSp], function(tmp) {
			frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
			names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]',
				unlist(strsplit(tmp,split="")))]
			return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
		}))

		locRead1= sam1[,"pos"]
		chrRead1= sam1[,"rName"]
		chrRead2= c()



		frqEventRead1 =lapply(sam1IndSp, function (tmp) {
			frqEventTmp= as.numeric(unlist(strsplit(sam1$cigar[tmp],
				split="[A-Z]")))
			names(frqEventTmp)= unlist(strsplit(sam1$cigar[tmp], 
				split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp]
					, split="")))]
			frqEventTmp= frqEventTmp[!is.na(match(names(frqEventTmp),
				c("N","D","M")))]

	# Finding start and end of the mapped pieces of the reads using cigar 
			eventTmp<- c(locRead1[tmp], locRead1[tmp]+
				cumsum(as.numeric(frqEventTmp)))
			begEventTmp<-eventTmp[-length(eventTmp)]
			endEventTmp<- begEventTmp+as.numeric(frqEventTmp)-1

			frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
				"M"))]
			return( list(readNo=rep(tmp,length(frqEventTmpM)), 
				chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
				frqEventTmpM, begEventTmpM, endEventTmpM) )

		})


		frqEventRead2=NA

		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]]))))
		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[2]]))))
		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]]))))
		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]]))))
		noRead2= c()
		chrRead2= c()
		begRead2= c()
		endRead2= c()

# analyzing reads
		noVec=c(noRead1,noRead2)
		chrVec=c(as.character(chrRead1),
			as.character(chrRead2))
		begVec=c(begRead1,begRead2)
		endVec=c(endRead1,endRead2)

		ref=reference[which(referenceIntronExon=="intron"),]
		readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
			IRanges::IRanges(start=begVec, end=endVec))

		refIntBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
			IRanges::IRanges(start=ref[,"begin"]-1, end=ref[,"begin"]-1))
		refIntEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
			IRanges::IRanges(start=ref[,"end"]+1, end=ref[,"end"]+1))

		hitsBegInt <- GenomicRanges::findOverlaps(readGRanges, 
			refIntBegGRanges, type= "end")
		hitsEndInt <- GenomicRanges::findOverlaps(readGRanges, 
			refIntEndGRanges, type= "start")
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GenomicRanges::GRanges( 
					seqnames=repeatsTableToFilter[,"chr"], 
					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
						end=repeatsTableToFilter[,"end"], 
						width=repeatsTableToFilter[,"end"]-
							repeatsTableToFilter[,"begin"]+1))
				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
					repeatGRanges, type= "any")
				filtBegInd= which(is.na(match(S4Vectors::queryHits(hitsBegInt),
					S4Vectors::queryHits(hitsFilter))))
				hitsBegInt=hitsBegInt[filtBegInd,]
				filtEndInd= which(is.na(match(S4Vectors::queryHits(hitsEndInt),
					S4Vectors::queryHits(hitsFilter))))
				hitsEndInt=hitsEndInt[filtEndInd,]
			}
			# Filter reads mapped completely to Either introns or exons

			refRes= rep(0,nrow(ref))
			if(length(c(S4Vectors::queryHits(hitsBegInt), 
				S4Vectors::queryHits(hitsEndInt)))>0){
				hitsApply=tapply(
					noVec[c(S4Vectors::queryHits(hitsBegInt), 
						S4Vectors::queryHits(hitsEndInt))], 
					c(S4Vectors::subjectHits(hitsBegInt), 
						S4Vectors::subjectHits(hitsEndInt)),
					function(tmp) return(length(unique(tmp))) )
				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
			}
			intSpanRes<- rep(0,nrow(reference))
			intSpanRes[which(referenceIntronExon=="intron")]<- refRes
	}

########!!!
########
	finalRes<-c()
	if(length(intRetRes)>0)
		finalRes<-c(finalRes, intRetRes)
	if(length(exExRes)>0){
		finalRes<-c(finalRes, exExRes)
	}
########
########!!!
	if(length(intSpanRes)>0){
		finalRes<-c(finalRes, intSpanRes)
	}
########!!!
########
	return(finalRes)
#End of defining function for single mapped reads that runs on parallel cores
} 
	

