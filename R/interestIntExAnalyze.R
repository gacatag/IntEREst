interestIntExAnalyse <- function(
	no, 
	reference,
	bamPrerocessRes,
	bamFile,
	yieldSize,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly,
	includePairedDuplicate,
	includeSingleReadDuplicate)
{

	bamPrerocessResRowNo <- no
	bpCount<- as.numeric(unlist(
		bamPrerocessRes$datReadNo[bamPrerocessResRowNo]))
	paired<- unlist(
		bamPrerocessRes$isPaired[bamPrerocessResRowNo])
	#Define functions
	subjectHits<-  S4Vectors::subjectHits
	queryHits<- S4Vectors::queryHits
	GRanges<- GenomicRanges::GRanges
	findOverlaps<- GenomicRanges::findOverlaps
	IRanges<- IRanges::IRanges

	if(logFile!="")
		cat( "InERESt:interestAnalyse: Analyzing read set ",
			bamPrerocessResRowNo, " out of ", nrow(bamPrerocessRes), 
			"\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse: Analyzing read set ",
		bamPrerocessResRowNo, " out of ", nrow(bamPrerocessRes), "\n")
	if(paired){
		bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, 
			asMates=TRUE )
		count=0
		open(bf)
# Read the correct block of read
		while(count<bpCount) {
			count=count+1; 
			scParam=Rsamtools::ScanBamParam(
				what=Rsamtools::scanBamWhat()[c(1,
					3,5,8,13,9, 10, 6, 4, 14, 15)], 
				flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				isDuplicate=includePairedDuplicate))
				readTmp=GenomicAlignments::readGAlignmentPairs(bf, 
				param=scParam)

			readSizeChk=2*length(readTmp)
	
		}
		close(bf)
	} else {
			bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, 
				asMates=TRUE)
		open(bf)
		readSizeChk=yieldSize 
		count=0
		while(count<bpCount) {
			count=count+1; 
			scParam=Rsamtools::ScanBamParam(
				what=Rsamtools::scanBamWhat()[c(1,
					3,5,8,13,9, 10, 6, 4,14,15)], 
				flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE, 
					isPaired=TRUE, 
					isDuplicate=includeSingleReadDuplicate))
				readTmp<- GenomicAlignments::readGAlignmentsList(bf, 
					param=scParam)
			lenTmp<- unlist(sapply(readTmp, length))
			readTmp<- readTmp[which(lenTmp<= maxNoMappedReads)]
		} 
		close(bf)
	}
#Extract information of the mappred pieces of the reads
	if(paired){
		sam1<- as.data.frame(GenomicAlignments::first(readTmp))[,
			c("qname", "rname", "strand", "pos", "qual", "cigar")]
		sam2<- as.data.frame(GenomicAlignments::second(readTmp))[,
			c("qname", "rname", "strand", "pos", "qual", "cigar")]
		colnames(sam1)<- c("qName","rName","strand","pos","qual",
			"cigar")
		colnames(sam2)<- c("qName","rName","strand","pos","qual",
			"cigar")


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

		noRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]])))
		chrRead1= unlist(sapply(frqEventRead1, function(tmp)
			as.character(unlist(tmp[[2]]))))
		begRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]])))
		endRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]])))
		noRead2= unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[1]])))
		chrRead2= unlist(sapply(frqEventRead2, function(tmp)
			as.character(unlist(tmp[[2]]))))
		begRead2= unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[4]])))
		endRead2= unlist(sapply(frqEventRead2, function(tmp)
			unlist(tmp[[5]])))
# if !paired	
	} else {
	
		sam1<- as.data.frame(readTmp)[, c("qname", "rname", "strand", 
			"pos", "qual", "cigar")]
		colnames(sam1)<- c("qName","rName","strand","pos","qual",
			"cigar")
		sam2= NA
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
		noRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[1]])))
		chrRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[2]])))
		begRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[4]])))
		endRead1= unlist(sapply(frqEventRead1, function(tmp)
			unlist(tmp[[5]])))
		noRead2= c()
		chrRead2= c()
		begRead2= c()
		endRead2= c()
#The end of else if			
	}
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
		readGRanges=GRanges( seqnames=chrVec, IRanges(start=begVec, 
			end=endVec))
		refGRanges= GRanges( seqnames=ref[,"chr"], 
			IRanges(start=ref[,"begin"], end=ref[,"end"], 
				width=ref[,"end"]-ref[,"begin"]+1))
		hits1 <- findOverlaps(readGRanges, refGRanges, type= "start")
		hits2 <- findOverlaps(readGRanges, refGRanges, type= "end")
		hits3 <- findOverlaps(readGRanges, refGRanges, type= "equal")
		hitsSubject<- c()
		hitsQuery<- c()
		# Filter reads mapped to repeats regions if requested
		if(length(repeatsTableToFilter)!=0){
			repeatGRanges= GRanges( 
				seqnames=repeatsTableToFilter[,"chr"], 
				IRanges(start=repeatsTableToFilter[,"begin"], 
					end=repeatsTableToFilter[,"end"], 
					width=repeatsTableToFilter[,"end"]-
						repeatsTableToFilter[,"begin"]+1))
			hitsFilter= findOverlaps(readGRanges, repeatGRanges, 
				type= "any")

			filtInd1= which(is.na(match(queryHits(hits1), 
				queryHits(hitsFilter))))
			filtInd2= which(is.na(match(queryHits(hits2),  
				queryHits(hitsFilter))))
			filtInd3= which(is.na(match(queryHits(hits3),  
				queryHits(hitsFilter))))
			hits1= hits1[filtInd1,]
			hits2= hits2[filtInd2,]
			hits3= hits3[filtInd3,]
		}
		# Include the reads that completely map to exons if user asked
		if(!junctionReadsOnly){
			hitsExtra= GenomicRanges::findOverlaps(readGRanges, 
				refGRanges, type= "within")
			hitsSubject=subjectHits(hitsExtra)
			hitsQuery=queryHits(hitsExtra)
		}
	
		hitsSubject=c(hitsSubject, subjectHits(hits1),
			subjectHits(hits2),subjectHits(hits3))
		hitsQuery=c(hitsQuery, queryHits(hits1), queryHits(hits2),
			queryHits(hits3))
		refRes= rep(0,nrow(ref))
		hitsApply= tapply(noVec[hitsQuery], hitsSubject, function(tmp) 
			return(length(unique(tmp))) )
		refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
		exExRes<-rep(0, nrow(reference))
		exExRes[which(referenceIntronExon=="exon")]<- refRes

	}
	methodNo=which(method=="IntRet")
	if(length(methodNo)>0){

		ref=reference[which(referenceIntronExon=="intron"),]
		readGRanges=GRanges( seqnames=chrVec, IRanges(start=begVec, 
			end=endVec))
		refGRanges= GRanges( seqnames=ref[,"chr"], 
			IRanges(start=ref[,"begin"], end=ref[,"end"], 
				width=ref[,"end"]-ref[,"begin"]+1))
		hits <- findOverlaps(readGRanges, refGRanges, type= "any")
		# Filter reads mapped to repeats regions if requested
		if(length(repeatsTableToFilter)!=0){
			repeatGRanges= GRanges( 
				seqnames=repeatsTableToFilter[,"chr"], 
				IRanges(start=repeatsTableToFilter[,"begin"],
					end=repeatsTableToFilter[,"end"], 
					width=repeatsTableToFilter[,"end"]-
						repeatsTableToFilter[,"begin"]+1))
			hitsFilter= findOverlaps(readGRanges, repeatGRanges, 
				type= "any")
			filtInd= which(is.na(match(queryHits(hits), 
				queryHits(hitsFilter))))
			hits=hits[filtInd,]
		}
		# Filter reads mapped completely to Either introns or exons
		if(junctionReadsOnly){
			hitsFilter= findOverlaps(readGRanges, refGRanges, 
				type= "within")
			filtInd=which(is.na(match(queryHits(hits), 
				queryHits(hitsFilter))))
			hits=hits[filtInd,]
		}
		refRes= rep(0,nrow(ref))
		hitsApply=tapply(noVec[queryHits(hits)], subjectHits(hits), 
			function(tmp) return(length(unique(tmp))) )
		refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
		intRetRes<- rep(0,nrow(reference))
		intRetRes[which(referenceIntronExon=="intron")]<- refRes

	}
	finalRes<-c()
	if(length(intRetRes)>0)
		finalRes<-c(finalRes, intRetRes)
	if(length(exExRes)>0){
		finalRes<-c(finalRes, exExRes)
	}
	return(finalRes)
#End of defining function that runs on all parallel cores 
} 
	
