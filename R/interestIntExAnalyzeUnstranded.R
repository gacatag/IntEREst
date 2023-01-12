interestIntExAnalysePairUnstranded <- function(
	readTmp,
	reference,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly,
	limitRanges,
	excludeFusionReads)
{
  qNam<- as.data.frame(GenomicAlignments::first(readTmp))[1,"qname"]
  if(logFile!="")
  	cat( "start ", 
  	   qNam, "\n",
  		file=logFile, append=appendLogFile)
  cat( "start ", 
       qNam, "\n")

  if(length(limitRanges)>0){
    r1Map<- GenomicRanges::findOverlaps(
      GenomicAlignments::first(readTmp), 
      limitRanges,
      type= "within")
    r2Map<- GenomicRanges::findOverlaps(
      GenomicAlignments::second(readTmp), 
      limitRanges,
      type= "within")
    filLimInd<- sort(unique(intersect(S4Vectors::queryHits(r1Map),
                                      S4Vectors::queryHits(r2Map))), 
                     decreasing=FALSE)
    filExtra<- c()
    if(excludeFusionReads){
      r1LimReads<- tapply(queryHits(r1Map),subjectHits(r1Map),  unique)
      r2LimReads<- tapply(queryHits(r2Map),subjectHits(r2Map),  unique)
      commonSubjects<- intersect(names(r1LimReads), names(r2LimReads))
      filExtraTmp<- lapply(commonSubjects, 
                           function(x) return (
                             intersect(r1LimReads[[x]], r2LimReads[[x]])
                             ))
      filExtra<-  unique(unlist(filExtraTmp))
      if(length(filExtra[which(!is.na(filExtra))])>0)
        filLimInd<- intersect(filLimInd, filExtra)
    }
    readTmp<- readTmp[filLimInd,]
  }
	sam1<- as.data.frame(GenomicAlignments::first(readTmp))[,
		c("qname", "rname", "strand", "pos", "qual", "cigar")]
	sam2<- as.data.frame(GenomicAlignments::second(readTmp))[,
		c("qname", "rname", "strand", "pos", "qual", "cigar")]
	colnames(sam1)<- c("qName","rName","strand","pos","qual",
		"cigar")
	colnames(sam2)<- c("qName","rName","strand","pos","qual",
		"cigar")
	# if(logFile!="")
	# 	cat( "Analyzing paired reads, ID: ", sam1[1,"qName"], "\n", 
	# 		file=logFile, append=appendLogFile)
	# cat( "Analyzing paired reads, ID: ", sam1[1,"qName"], "\n")
######
######!!!
	exExRes<- c()
	intRetRes<- c()
	intSpanRes<- c()
	exSkipRes<- c()
  if(length(readTmp)>0)	{
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
  			
  			if(logFile!="")
  			  cat( "Part1 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
  			
  			hits1 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "start")
  			
  			if(logFile!="")
  			  cat( "Part2 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits2 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "end")
  			
  			if(logFile!="")
  			  cat( "Part3 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits3 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "equal")
  			
  			if(logFile!="")
  			  cat( "Part4 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
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
  			
  			if(logFile!="")
  			  cat( "Part5 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
  			hitsQuery=c(hitsQuery, S4Vectors::queryHits(hits1), 
  				S4Vectors::queryHits(hits2), S4Vectors::queryHits(hits3))
  			
  			if(logFile!="")
  			  cat( "Part6 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			refRes= rep(0,nrow(ref))
  			if(length(hitsQuery)>0){
  				hitsApply= tapply(noVec[hitsQuery], hitsSubject, function(tmp)
  					return(length(unique(tmp))) )
  				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  			}
  			if(logFile!="")
  			  cat( "Part7 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			exExRes<-rep(0, nrow(reference))
  			exExRes[which(referenceIntronExon=="exon")]<- refRes
  			rm("refRes")
  	
  		}	
  		methodNo=which(method=="IntRet")
  		if((length(methodNo)>0) & (! junctionReadsOnly)){
  
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
  			refRes= rep(0,nrow(ref))
  			if(length(S4Vectors::queryHits(hits))>0){
  				hitsApply=tapply(noVec[S4Vectors::queryHits(hits)], 
  					S4Vectors::subjectHits(hits), 
  					function(tmp) return(length(unique(tmp))) )
  				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  			}
  			intRetRes<- rep(0,nrow(reference))
  			intRetRes[which(referenceIntronExon=="intron")]<- refRes
  			rm("refRes")
  	
  		} else if ((length(methodNo)>0) & (junctionReadsOnly)) {
  			ref=reference[which(referenceIntronExon=="intron"),]
  			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  				IRanges::IRanges(start=begVec, end=endVec))
  			refEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"end"], end=ref[,"end"]))
  			refBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"begin"]))
  			hitsEnd <- GenomicRanges::findOverlaps(refEndGRanges, readGRanges,
  				type= "within")
  			hitsBeg <- GenomicRanges::findOverlaps(refBegGRanges, readGRanges,
  				type= "within")
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
  				filtEndInd= which(is.na(match(S4Vectors::subjectHits(hitsEnd), 
  					S4Vectors::queryHits(hitsFilter))))
  				filtBegInd= which(is.na(match(S4Vectors::subjectHits(hitsBeg), 
  					S4Vectors::queryHits(hitsFilter))))
  				hitsEnd=hitsEnd[filtEndInd,]
  				hitsBeg=hitsBeg[filtBegInd,]
  			}
  
  			refRes= rep(0,nrow(ref))
  			refResBeg= refRes
  			refResEnd= refRes
  			if( (length(S4Vectors::subjectHits(hitsBeg))>0) | 
  				(length(S4Vectors::subjectHits(hitsEnd))>0) ){
  				hitsEndApply=tapply(noVec[S4Vectors::subjectHits(hitsEnd)], 
  					S4Vectors::queryHits(hitsEnd), 
  					function(tmp) return(length(unique(tmp))) )
  				hitsBegApply=tapply(noVec[S4Vectors::subjectHits(hitsBeg)], 
  					S4Vectors::queryHits(hitsBeg), 
  					function(tmp) return(length(unique(tmp))) )
  				refResBeg[as.numeric(names(hitsBegApply))]= 
  					as.numeric(hitsBegApply)
  				refResEnd[as.numeric(names(hitsEndApply))]= 
  					as.numeric(hitsEndApply)
  				refRes=refResBeg+refResEnd
  				rm("refResBeg","refResEnd")
  			}
  			intRetRes<- rep(0,nrow(reference))
  			intRetRes[which(referenceIntronExon=="intron")]<- refRes
  			rm("refRes")
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
  			rm("refRes")
  	}
  
  	####
  	### NEW ExSkip
  	methodNo=which(method=="ExSkip")
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
  	  
  	  locRead1= sam1[,"pos"]
  	  chrRead1= sam1[,"rName"]
  	  locRead2= sam2[,"pos"]
  	  chrRead2= sam2[,"rName"]
  	  
  	  
  	  
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
  	    frqEventTmpN=frqEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    begEventTmpN=begEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    endEventTmpN=endEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    return( list(readNo=rep(tmp,length(frqEventTmpN)), 
  	                 chr=rep(chrRead1[tmp],length(frqEventTmpN)), 
  	                 frqEventTmpN, begEventTmpN, endEventTmpN,
  	                 frqEventTmpM, begEventTmpM, endEventTmpM) )
  	    
  	  })
  	  
  	  
  	  frqEventRead2= lapply(sam2IndSp, function (tmp) {
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
  	    frqEventTmpN=frqEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    begEventTmpN=begEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    endEventTmpN=endEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    return( list(readNo=rep(tmp,length(frqEventTmpN)), 
  	                 chr=rep(chrRead1[tmp],length(frqEventTmpN)), 
  	                 frqEventTmpN, begEventTmpN, endEventTmpN,
  	                 frqEventTmpM, begEventTmpM, endEventTmpM,
  	                 readNoM=rep(tmp,length(frqEventTmpM))) )
  	    
  	  })
  	  
  	  noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[1]]))))
  	  chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[2]]))))
  	  begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[4]]))))
  	  endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[5]]))))
  
  	  noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[1]]))))
  	  chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[2]]))))
  	  begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[4]]))))
  	  endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[5]]))))
  	  
  	  # analyzing reads
  	  noVec=c(noRead1,noRead2)
  	  chrVec=c(as.character(chrRead1),
  	           as.character(chrRead2))
  	  begVec=c(begRead1,begRead2)
  	  endVec=c(endRead1,endRead2)
  	  
  	  ref=reference[which(referenceIntronExon=="exon"),]
  	  readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  	                                      IRanges::IRanges(start=begVec, end=endVec))
  	  
  	  refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  	                                      IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"]))
  	  
  	  hitsAllEx <- GenomicRanges::findOverlaps(refGRanges, 
  	                                           readGRanges, type= "within")
  	  
  	  # Filter reads mapped to repeats regions if requested
  	  if(length(repeatsTableToFilter)!=0){
  	    noRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[9]]))))
  	    chrRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[2]]))))
  	    begRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[7]]))))
  	    endRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[8]]))))
  
  	    noRead2M= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	      unlist(tmp[[9]]))))
  	    chrRead2M= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	      unlist(tmp[[2]]))))
  	    begRead2M= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	      unlist(tmp[[7]]))))
  	    endRead2M= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	      unlist(tmp[[8]]))))
  	    
  	    # analyzing reads
  	    noVecM=c(noRead1M,noRead2M)
  	    chrVecM=c(as.character(chrRead1M),
  	              as.character(chrRead2M))
  	    begVecM=c(begRead1M,begRead2M)
  	    endVecM=c(endRead1M,endRead2M)
  	    
  	    readGRangesM=GenomicRanges::GRanges( seqnames=chrVecM, 
  	                                         IRanges::IRanges(start=begVecM, end=endVecM))
  	    repeatGRanges= GenomicRanges::GRanges( 
  	      seqnames=repeatsTableToFilter[,"chr"], 
  	      IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  	                       end=repeatsTableToFilter[,"end"], 
  	                       width=repeatsTableToFilter[,"end"]-
  	                         repeatsTableToFilter[,"begin"]+1))
  	    hitsFilter= GenomicRanges::findOverlaps(readGRangesM, 
  	                                            repeatGRanges, type= "any")
  	    # filthitsAllExInd= which(is.na(match(S4Vectors::subjectHits(hitsBegInt),
  	    #                               S4Vectors::queryHits(hitsFilter))))
  	    noVecFil<- unique(noVecM[queryHits(hitsFilter)])
  	    filthitsAllExInd<- which(!subjectHits(hitsAllEx) %in% 
  	                               subjectHits(hitsAllEx)[
  	                                 which(noVec[subjectHits(hitsAllEx)]
  	                                       %in% noVecFil)])
  	    hitsAllEx=hitsAllEx[filthitsAllExInd,]
  	    
  	  }
  	  # Filter reads mapped completely to Either introns or exons
  	  
  	  refRes= rep(0,nrow(ref))
  	  if(length(S4Vectors::subjectHits(hitsAllEx))>0){
  	    hitsApply=tapply(
  	      noVec[S4Vectors::subjectHits(hitsAllEx)], 
  	      S4Vectors::queryHits(hitsAllEx),
  	      function(tmp) return(length(unique(tmp))) )
  	    refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  	  }
  	  exSkipRes<- rep(0,nrow(reference))
  	  exSkipRes[which(referenceIntronExon=="exon")]<- refRes
  	  rm("refRes")
  	}
  }
########!!!
########	
	
########!!!
########
	finalRes<-c()
	
	if(length(exExRes)==0)
	  exExRes<- rep(0,nrow(reference))
	
	if(length(intRetRes)==0)
	  intRetRes<- rep(0,nrow(reference))
	  
	if(length(intSpanRes)==0)
	  intSpanRes<- rep(0,nrow(reference))

	if(length(exSkipRes)==0)
	  exSkipRes<- rep(0,nrow(reference))

	#	if(length(exExRes)>0){
	if("ExEx" %in% method)
	  finalRes<-c(finalRes, exExRes)
	rm("exExRes")
	
#	if(length(intRetRes)>0)
	if("IntRet" %in% method)
		finalRes<-c(finalRes, intRetRes)
	rm("intRetRes")

########
########!!!
#	if(length(intSpanRes)>0)
	if("IntSpan" %in% method)
		finalRes<-c(finalRes, intSpanRes)
	rm("intSpanRes")
	
#	if(length(exSkipRes)>0)
	if("ExSkip" %in% method)
	  finalRes<-c(finalRes, exSkipRes)
	rm("exSkipRes")
########!!!
########
	if(logFile!="")
	  cat( "end ", 
	       qNam, " ", length(finalRes), "\n",
	       file=logFile, append=appendLogFile)
	cat( "end ", 
	     qNam, " ", length(finalRes), "\n")	
	# if(is.null(finalRes))
	#   finalRes<- rep(0, 
	#     length(which(method%in%c("ExEx", "IntRet", "IntSpan", "ExSkip"))))
	return(finalRes)
	rm("finalRes")
#End of defining function for paired-mapped reads that runs on parallel cores 
} 
	

# For single reads
interestIntExAnalyseSingleUnstranded <- function(
	readTmp,
	reference,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly,
	limitRanges,
	excludeFusionReads)
{
  qNam<- as.data.frame(readTmp)[1, "qname"]
  if(logFile!="")
    cat( "Singlestart ", 
         qNam, "\n",
         file=logFile, append=appendLogFile)
  cat( "Singlestart ", 
       qNam, "\n")
  
# Filtering readTmp based on number of mapping targets for each read
	lenTmp<- unlist(sapply(readTmp, length))
	readTmpFil<- readTmp[which(lenTmp<= maxNoMappedReads)]

	if(length(limitRanges)>0 & length(readTmpFil)>0){
	  r1Map<- GenomicRanges::findOverlaps(
	    readTmpFil, 
	    limitRanges,
	    type= "within")
	  
	  filLimInd<- sort(unique(S4Vectors::queryHits(r1Map)), decreasing=FALSE)
	  readTmpFil<- readTmpFil[filLimInd,]
	}
	
#Extract information of the mapped pieces of the reads
	sam1<- as.data.frame(readTmpFil)[, c("qname", "rname", "strand", 
		"pos", "qual", "cigar")]
	colnames(sam1)<- c("qName","rName","strand","pos","qual",
		"cigar")
	sam2= NA

	# if(logFile!="")
	# 	cat( "Analyzing single reads, ID: ", sam1[1,"qName"], "\n", 
	# 		file=logFile, append=appendLogFile)
	# cat( "Analyzing single reads, ID: ", sam1[1,"qName"], "\n")

######
######!!!
	exExRes<- c()
	intRetRes<- c()
	intSpanRes<- c()
	exSkipRes<- c()
  if(length(readTmp)>0){
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
  		methodNo= which(method=="ExEx")
  		if(length(methodNo)>0){
  			ref= reference[which(referenceIntronExon=="exon"),]
  			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  				IRanges::IRanges(start=begVec, end=endVec))
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
  					width=ref[,"end"]-ref[,"begin"]+1))
  			hits1 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "start")
  			if(logFile!="")
  			  cat( "Part8 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			hits2 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "end")
  			if(logFile!="")
  			  cat( "Part9 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			hits3 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "equal")
  			if(logFile!="")
  			  cat( "Part10 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
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
  			if(logFile!="")
  			  cat( "Part11 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			hitsQuery=c(hitsQuery, S4Vectors::queryHits(hits1), 
  				S4Vectors::queryHits(hits2), S4Vectors::queryHits(hits3))
  			refRes= rep(0,nrow(ref))
  			if(length(hitsQuery)>0){
  				hitsApply= tapply(noVec[hitsQuery], hitsSubject, function(tmp)
  					return(length(unique(tmp))) )
  				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  			}
  			if(logFile!="")
  			  cat( "Part12 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			exExRes<-rep(0, nrow(reference))
  			exExRes[which(referenceIntronExon=="exon")]<- refRes

  			rm("refRes")
  	
  		}
  		methodNo=which(method=="IntRet")
  		if((length(methodNo)>0) & (! junctionReadsOnly)){
  	
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
  			refRes= rep(0,nrow(ref))
  			if(length(S4Vectors::queryHits(hits))>0){
  				hitsApply=tapply(noVec[S4Vectors::queryHits(hits)], 
  					S4Vectors::subjectHits(hits), 
  					function(tmp) return(length(unique(tmp))) )
  				refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  			}
  			intRetRes<- rep(0,nrow(reference))
  			intRetRes[which(referenceIntronExon=="intron")]<- refRes
  			
  			rm("refRes")
  	
  		} else if ((length(methodNo)>0) & (junctionReadsOnly)) {
  			ref=reference[which(referenceIntronExon=="intron"),]
  			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  				IRanges::IRanges(start=begVec, end=endVec))
  			refEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"end"], end=ref[,"end"]))
  			refBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"begin"]))
  			hitsEnd <- GenomicRanges::findOverlaps(refEndGRanges, readGRanges,
  				type= "within")
  			hitsBeg <- GenomicRanges::findOverlaps(refBegGRanges, readGRanges,
  				type= "within")
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
  				filtEndInd= which(is.na(match(S4Vectors::subjectHits(hitsEnd), 
  					S4Vectors::queryHits(hitsFilter))))
  				filtBegInd= which(is.na(match(S4Vectors::subjectHits(hitsBeg), 
  					S4Vectors::queryHits(hitsFilter))))
  				hitsEnd=hitsEnd[filtEndInd,]
  				hitsBeg=hitsBeg[filtBegInd,]
  			}
  
  			refRes= rep(0,nrow(ref))
  			refResBeg= refRes
  			refResEnd= refRes
  			if( (length(S4Vectors::subjectHits(hitsBeg))>0) | 
  				(length(S4Vectors::subjectHits(hitsEnd))>0) ){
  				hitsEndApply=tapply(noVec[S4Vectors::subjectHits(hitsEnd)], 
  					S4Vectors::queryHits(hitsEnd), 
  					function(tmp) return(length(unique(tmp))) )
  				hitsBegApply=tapply(noVec[S4Vectors::subjectHits(hitsBeg)], 
  					S4Vectors::queryHits(hitsBeg), 
  					function(tmp) return(length(unique(tmp))) )
  				refResBeg[as.numeric(names(hitsBegApply))]= 
  					as.numeric(hitsBegApply)
  				refResEnd[as.numeric(names(hitsEndApply))]= 
  					as.numeric(hitsEndApply)
  				refRes=refResBeg+refResEnd
  				rm("refResBeg","refResEnd")
  			}
  			intRetRes<- rep(0,nrow(reference))
  			intRetRes[which(referenceIntronExon=="intron")]<- refRes
  			
  			rm("refRes")
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
  			
  			rm("refRes")
  	}
  
  ####
  ### NEW ExSkip
  	methodNo=which(method=="ExSkip")
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
  	    frqEventTmpN=frqEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    begEventTmpN=begEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    endEventTmpN=endEventTmp[!is.na(match(names(frqEventTmp),
  	                                          "N"))]
  	    return( list(readNo=rep(tmp,length(frqEventTmpN)), 
  	                 chr=rep(chrRead1[tmp],length(frqEventTmpN)), 
  	                 frqEventTmpN, begEventTmpN, endEventTmpN,
  	                 frqEventTmpM, begEventTmpM, endEventTmpM,
  	                 readNoM=rep(tmp,length(frqEventTmpM))) )
  	    
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
  	  
  	  ref=reference[which(referenceIntronExon=="exon"),]
  	  readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  	                                      IRanges::IRanges(start=begVec, end=endVec))
  	  
  	  refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  	                                            IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"]))
  
  	  hitsAllEx <- GenomicRanges::findOverlaps(refGRanges, 
  	                                           readGRanges, type= "within")
  
  	  # Filter reads mapped to repeats regions if requested
  	  if(length(repeatsTableToFilter)!=0){
  	    noRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[9]]))))
  	    chrRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[2]]))))
  	    begRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[7]]))))
  	    endRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[8]]))))
  	    noRead2M= c()
  	    chrRead2M= c()
  	    begRead2M= c()
  	    endRead2M= c()
  	    
  	    # analyzing reads
  	    noVecM=c(noRead1M,noRead2M)
  	    chrVecM=c(as.character(chrRead1M),
  	             as.character(chrRead2M))
  	    begVecM=c(begRead1M,begRead2M)
  	    endVecM=c(endRead1M,endRead2M)
  	    
  	    readGRangesM=GenomicRanges::GRanges( seqnames=chrVecM, 
  	                                        IRanges::IRanges(start=begVecM, end=endVecM))
  	    repeatGRanges= GenomicRanges::GRanges( 
  	      seqnames=repeatsTableToFilter[,"chr"], 
  	      IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  	                       end=repeatsTableToFilter[,"end"], 
  	                       width=repeatsTableToFilter[,"end"]-
  	                         repeatsTableToFilter[,"begin"]+1))
  	    hitsFilter= GenomicRanges::findOverlaps(readGRangesM, 
  	                                            repeatGRanges, type= "any")
  	    # filthitsAllExInd= which(is.na(match(S4Vectors::subjectHits(hitsBegInt),
  	    #                               S4Vectors::queryHits(hitsFilter))))
  	    noVecFil<- unique(noVecM[queryHits(hitsFilter)])
  	    filthitsAllExInd<- which(!subjectHits(hitsAllEx) %in% 
  	                               subjectHits(hitsAllEx)[
  	                                 which(noVec[subjectHits(hitsAllEx)]
  	                                       %in% noVecFil)])
  	    hitsAllEx=hitsAllEx[filthitsAllExInd,]
  
  	  }
  	  # Filter reads mapped completely to Either introns or exons
  	  
  	  refRes= rep(0,nrow(ref))
  	  if(length(S4Vectors::subjectHits(hitsAllEx))>0){
  	    hitsApply=tapply(
  	      noVec[S4Vectors::subjectHits(hitsAllEx)], 
  	      S4Vectors::queryHits(hitsAllEx),
  	      function(tmp) return(length(unique(tmp))) )
  	    refRes[as.numeric(names(hitsApply))]= as.vector(hitsApply)
  	  }
  	  exSkipRes<- rep(0,nrow(reference))
  	  exSkipRes[which(referenceIntronExon=="exon")]<- refRes
  	  
  	  rm("refRes")
  	}
  }
########!!!
########
	finalRes<-c()
	
	if(length(exExRes)==0)
	  exExRes<- rep(0,nrow(reference))
	
	if(length(intRetRes)==0)
	  intRetRes<- rep(0,nrow(reference))
	
	if(length(intSpanRes)==0)
	  intSpanRes<- rep(0,nrow(reference))
	
	if(length(exSkipRes)==0)
	  exSkipRes<- rep(0,nrow(reference))

	#	if(length(exExRes)>0){
	if("ExEx" %in% method)
	  finalRes<-c(finalRes, exExRes)
	rm("exExRes")	
		
	#	if(length(intRetRes)>0)
	if("IntRet" %in% method)
	  finalRes<-c(finalRes, intRetRes)
	rm("intRetRes")
	
	########
	########!!!
	#	if(length(intSpanRes)>0)
	if("IntSpan" %in% method)
	  finalRes<-c(finalRes, intSpanRes)
	rm("intSpanRes")
	
	#	if(length(exSkipRes)>0)
	if("ExSkip" %in% method)
	  finalRes<-c(finalRes, exSkipRes)
	rm("exSkipRes")
########!!!
########
	if(logFile!="")
	  cat( "Singleend ", 
	       qNam, " ", length(finalRes), "\n",
	       file=logFile, append=appendLogFile)
	cat( "Singleend ", 
	     qNam, " ", length(finalRes), "\n")
	# if(is.null(finalRes) | (length(finalRes)==0))
	#   finalRes<- rep(0, 
	#     length(which(method%in%c("ExEx", "IntRet", "IntSpan", "ExSkip"))))
	return(finalRes)
	rm("finalRes")
#End of defining function for single mapped reads that runs on parallel cores
} 
	

