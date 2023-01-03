interestIntExAnalysePairStranded <- function(
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
	strandSpecific,
	excludeFusionReads){
  
  revStr<- c("-", "+", "*")
  names(revStr)<- c("+", "-", "*") 
  qNam<- as.data.frame(GenomicAlignments::first(readTmp))[1,"qname"]
  if(logFile!="")
  	cat( "start ", 
  	   qNam, "\n",
  		file=logFile, append=appendLogFile)
  cat( "start ", 
       qNam, "\n")
  
  if(length(limitRanges)>0){
    
    readTmpFirst<- GenomicAlignments::first(readTmp)
    readTmpSecond<- GenomicAlignments::second(readTmp)
    if(strandSpecific=="stranded"){
      GenomicRanges::strand(readTmpFirst)<- 
        as.character(GenomicRanges::strand(readTmpFirst))
      GenomicRanges::strand(readTmpSecond)<- 
        revStr[as.character(GenomicRanges::strand(readTmpSecond))]
    } else if (strandSpecific=="reverse"){
      GenomicRanges::strand(readTmpSecond)<- 
        as.character(GenomicRanges::strand(readTmpSecond))
      GenomicRanges::strand(readTmpFirst)<- 
        revStr[as.character(GenomicRanges::strand(readTmpFirst))]
    }
    
    r1Map<- GenomicRanges::findOverlaps(
      readTmpFirst, 
      limitRanges,
      type= "within",
      ignore.strand= FALSE)
    r2Map<- GenomicRanges::findOverlaps(
      readTmpSecond, 
      limitRanges,
      type= "within",
      ignore.strand= FALSE)
    
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
  
  filStrandInd<-which ( xor(
    (as.data.frame(GenomicAlignments::first(readTmp))[,"strand"]=="+") , 
    (as.data.frame(GenomicAlignments::second(readTmp))[,"strand"]=="+") ))
  readTmp<- readTmp[filStrandInd,]
  
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
  		strRead1=sam1[,"strand"]
  		chrRead2=sam2[,"rName"]
  		locRead2=sam2[,"pos"]
  		strRead2=sam2[,"strand"]
  		if(strandSpecific=="stranded"){
  		  strAdj1<- as.character(strRead1)
  		  strAdj2<- as.character(revStr[strRead2])
  		} else if (strandSpecific=="reverse"){
  		  strAdj2<- as.character(strRead2)
  		  strAdj1<- as.character(revStr[strRead1])
  		}
  		
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
  				frqEventTmpM, begEventTmpM, endEventTmpM, 
  				str=rep(strRead1[tmp],length(frqEventTmpM))) )
  	
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
  				frqEventTmpM, begEventTmpM, endEventTmpM, 
  				str=rep(strRead2[tmp],length(frqEventTmpM))))
  		})
  
  		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[1]]))))
  		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			as.character(unlist(tmp[[2]])))))
  		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[4]]))))
  		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[5]]))))
  		strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[1]]))))
  		chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			as.character(unlist(tmp[[2]])))))
  		begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[4]]))))
  		endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[5]]))))
  		strRead2=as.vector(unlist(sapply(frqEventRead2, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		
  		if(strandSpecific=="stranded"){
  		  strAdj1<- as.character(strRead1)
  		  strAdj2<- as.character(revStr[strRead2])
  		} else if (strandSpecific=="reverse"){
  		  strAdj2<- as.character(strRead2)
  		  strAdj1<- as.character(revStr[strRead1])
  		}
  		
  		noVec=c(noRead1,noRead2)
  		chrVec=c(as.character(chrRead1),
  			as.character(chrRead2))
  		begVec=c(begRead1,begRead2)
  		endVec=c(endRead1,endRead2)
  		strVec=c(as.character(strAdj1),
  		         as.character(strAdj2))
  
  		methodNo=which(method=="ExEx")
  		if(length(methodNo)>0){
  			ref=reference[which(referenceIntronExon=="exon"),]
  			
  			#if(strandSpecific%in%c("stranded","reverse")){
  			# FOR READ1
  			read1GRanges=GenomicRanges::GRanges( seqnames=chrRead1, 
  				IRanges::IRanges(start=begRead1, end=endRead1), strand=strAdj1)
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
  					width=ref[,"end"]-ref[,"begin"]+1), strand=ref[,"strand"])
  			
  			if(logFile!="")
  			  cat( "Part1 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
  			
  			hits1_1 <- GenomicRanges::findOverlaps(read1GRanges, refGRanges, 
  				type= "start", ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part2 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits1_2 <- GenomicRanges::findOverlaps(read1GRanges, refGRanges, 
  				type= "end", ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part3 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits1_3 <- GenomicRanges::findOverlaps(read1GRanges, refGRanges, 
  				type= "equal", ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part4 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			# FOR READ2
  			
  			read2GRanges=GenomicRanges::GRanges( seqnames=chrRead2, 
  			  IRanges::IRanges(start=begRead2, 
  			    end=endRead2), 
  			  strand=strAdj2)
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			  IRanges::IRanges(start=ref[,"begin"], 
  			    end=ref[,"end"], 
  			      width=ref[,"end"]-ref[,"begin"]+1), 
  			  strand=ref[,"strand"])
  			
  			if(logFile!="")
  			  cat( "Part1 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
  			
  			hits2_1 <- GenomicRanges::findOverlaps(read2GRanges, refGRanges, 
  			  type= "start", ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part2 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits2_2 <- GenomicRanges::findOverlaps(read2GRanges, refGRanges, 
  			                                       type= "end", 
  			                                       ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part3 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			hits2_3 <- GenomicRanges::findOverlaps(read2GRanges, refGRanges, 
  			                                       type= "equal", 
  			                                       ignore.strand=FALSE)
  			
  			if(logFile!="")
  			  cat( "Part4 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile) 			
  			hitsSubject<- c()
  			hitsQuery<- c()
  			hits1Subject=c()
  			hits1Query=c()
  			hits2Subject=c()
  			hits2Query=c()
  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"], 
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter= GenomicRanges::findOverlaps(read1GRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
  	
  				filtInd1_1= which(is.na(match(S4Vectors::queryHits(hits1_1), 
  					S4Vectors::queryHits(hitsFilter))))
  				filtInd1_2= which(is.na(match(S4Vectors::queryHits(hits1_2),  
  					S4Vectors::queryHits(hitsFilter))))
  				filtInd1_3= which(is.na(match(S4Vectors::queryHits(hits1_3),  
  					S4Vectors::queryHits(hitsFilter))))
  				
  				hitsFilter= GenomicRanges::findOverlaps(read2GRanges, 
  				                                        repeatGRanges, type= "any", ignore.strand=TRUE)
  				
  				filtInd2_1= which(is.na(match(S4Vectors::queryHits(hits2_1), 
  				                              S4Vectors::queryHits(hitsFilter))))
  				filtInd2_2= which(is.na(match(S4Vectors::queryHits(hits2_2),  
  				                              S4Vectors::queryHits(hitsFilter))))
  				filtInd2_3= which(is.na(match(S4Vectors::queryHits(hits2_3),  
  				                              S4Vectors::queryHits(hitsFilter))))
  				
  				hits1_1= hits1_1[filtInd1_1,]
  				hits1_2= hits1_2[filtInd1_2,]
  				hits1_3= hits1_3[filtInd1_3,]
  				hits2_1= hits2_1[filtInd2_1,]
  				hits2_2= hits2_2[filtInd2_2,]
  				hits2_3= hits2_3[filtInd2_3,]
  			}
  			# Include the reads that completely map to exons if user asked
  			if(!junctionReadsOnly){
  				hits1Extra= GenomicRanges::findOverlaps(read1GRanges, 
  					refGRanges, type= "within", ignore.strand=FALSE)
  				hits2Extra= GenomicRanges::findOverlaps(read2GRanges, 
  				  refGRanges, type= "within", ignore.strand=FALSE)

  				hits1Subject=S4Vectors::subjectHits(hits1Extra)
  				hits1Query=S4Vectors::queryHits(hits1Extra)
  				hits2Subject=S4Vectors::subjectHits(hits2Extra)
  				hits2Query=S4Vectors::queryHits(hits2Extra)
  			}
  		
  			hits1Subject=c(hits1Subject, S4Vectors::subjectHits(hits1_1),
  				S4Vectors::subjectHits(hits1_2),S4Vectors::subjectHits(hits1_3))
  			hits2Subject=c(hits2Subject, S4Vectors::subjectHits(hits2_1),
  			   S4Vectors::subjectHits(hits2_2),S4Vectors::subjectHits(hits2_3))
  			
  			if(logFile!="")
  			  cat( "Part5 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			
  			hits1Query=c(hits1Query, S4Vectors::queryHits(hits1_1), 
  			  S4Vectors::queryHits(hits1_2), S4Vectors::queryHits(hits1_3))
  			hits2Query=c(hits2Query, S4Vectors::queryHits(hits2_1), 
  			  S4Vectors::queryHits(hits2_2), S4Vectors::queryHits(hits2_3))
  			
  			if(logFile!="")
  			  cat( "Part6 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			
  			refRes= rep(0,nrow(ref))
  			if(length(hits1Query)>0){
  				hitsApply= tapply(c(noRead1[hits1Query], noRead2[hits2Query]), 
  				  c(hits1Subject, hits2Subject), 
  				  function(tmp) return(length(unique(tmp))) )
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
  			read1GRanges=GenomicRanges::GRanges( seqnames=chrRead1, 
  			  IRanges::IRanges(start=begRead1, end=endRead1), strand=strAdj1)
  			read2GRanges=GenomicRanges::GRanges( seqnames=chrRead2, 
  			  IRanges::IRanges(start=begRead2, end=endRead2), strand=strAdj2)
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			  IRanges::IRanges(start=ref[,"begin"], 
  			    end=ref[,"end"], 
  			    width=ref[,"end"]-ref[,"begin"]+1), 
  			  strand=ref[,"strand"])
  			
  			hits1 <- GenomicRanges::findOverlaps(read1GRanges, refGRanges, 
  				type= "any", ignore.strand=FALSE)
  			hits2 <- GenomicRanges::findOverlaps(read2GRanges, refGRanges, 
  			  type= "any", ignore.strand=FALSE)

  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter1= GenomicRanges::findOverlaps(read1GRanges, 
  					repeatGRanges, type= "any",ignore.strand=TRUE)
  				hitsFilter2= GenomicRanges::findOverlaps(read2GRanges, 
  				  repeatGRanges, type= "any", ignore.strand=TRUE)
# Explanation: If both reads have regions overlapping repeat regions ignore them!
  				filNoHits<- 
  				  unique(intersect(noRead1[queryHits(hitsFilter1)], 
  				            noRead2[queryHits(hitsFilter2)]))
  				  				
  				filtInd1= 
  				  which( !(noRead1[S4Vectors::queryHits(hits1)] %in% filNoHits))
  				filtInd2= 
  				  which( !(noRead2[S4Vectors::queryHits(hits2)] %in% filNoHits))
  				hits1=hits1[filtInd1,]
  				hits2=hits2[filtInd2,]
  			}
  			refRes= rep(0,nrow(ref))
  			if(length(S4Vectors::queryHits(hits1))>0 | length(S4Vectors::queryHits(hits2))>0){
  				hitsApply=tapply(
  				  c(noRead1[S4Vectors::queryHits(hits1)],
  				    noRead2[S4Vectors::queryHits(hits2)]),
  					c(S4Vectors::subjectHits(hits1), S4Vectors::subjectHits(hits2)),
  					function(tmp) return(length(unique(tmp))) )
  				
  				refRes[as.numeric(unique(names(hitsApply)))]= 
  				  as.vector(hitsApply)
  			}
  			intRetRes<- rep(0,nrow(reference))
  			intRetRes[which(referenceIntronExon=="intron")]<- refRes
  			rm("refRes")

  		} else if ((length(methodNo)>0) & (junctionReadsOnly)) {
  		  
  		  
  			ref=reference[which(referenceIntronExon=="intron"),]
  			refEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			  IRanges::IRanges(start=ref[,"end"], end=ref[,"end"]), 
  			  strand=ref[,"strand"])
  			refBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			  IRanges::IRanges(start=ref[,"begin"], end=ref[,"begin"]),
  			  strand=ref[,"strand"])
  			
  			read1GRanges=GenomicRanges::GRanges( seqnames=chrRead1, 
  			  IRanges::IRanges(start=begRead1, end=endRead1), strand=strAdj1)
  			hitsEnd1 <- GenomicRanges::findOverlaps(refEndGRanges, read1GRanges,
  			  type= "within", ignore.strand=FALSE)
  			hitsBeg1 <- GenomicRanges::findOverlaps(refBegGRanges, read1GRanges,
  			  type= "within", ignore.strand=FALSE)

  			read2GRanges=GenomicRanges::GRanges( seqnames=chrRead2, 
  			  IRanges::IRanges(start=begRead2, end=endRead2), strand=strAdj2)
  			hitsEnd2 <- GenomicRanges::findOverlaps(refEndGRanges, read2GRanges,
  			  type= "within", ignore.strand=FALSE)
  			hitsBeg2 <- GenomicRanges::findOverlaps(refBegGRanges, read2GRanges,
  			  type= "within", ignore.strand=FALSE)

  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter1= GenomicRanges::findOverlaps(read1GRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
  				hitsFilter2= GenomicRanges::findOverlaps(read2GRanges, 
  				  repeatGRanges, type= "any", ignore.strand=TRUE)
  				
  				noFilQuHits<- 
  				  intersect(noRead1[queryHits(hitsFilter1)], 
  				             noRead2[queryHits(hitsFilter2)])
  				
  				filtEndInd1= which(! noRead1[S4Vectors::subjectHits(hitsEnd1)] %in% 
  				                     noFilQuHits)
  				filtBegInd1= which(! noRead1[S4Vectors::subjectHits(hitsBeg1)] %in% 
  				                     noFilQuHits)
  				
  				filtEndInd2= which(! noRead2[S4Vectors::subjectHits(hitsEnd2)] %in% 
  				                     noFilQuHits)
  				filtBegInd2= which(! noRead2[S4Vectors::subjectHits(hitsBeg2)] %in% 
  				                     noFilQuHits)
  				
  				
  				hitsEnd1=hitsEnd1[filtEndInd1,]
  				hitsBeg1=hitsBeg1[filtBegInd1,]
  				hitsEnd2=hitsEnd2[filtEndInd2,]
  				hitsBeg2=hitsBeg2[filtBegInd2,]
  			}
  			
  			  
  			refRes= rep(0,nrow(ref))
  			refResBeg= refRes
  			refResEnd= refRes
  			if( (length(S4Vectors::subjectHits(hitsBeg1))>0) | 
  				(length(S4Vectors::subjectHits(hitsEnd1))>0) | 
  				(length(S4Vectors::subjectHits(hitsBeg2))>0) | 
  				(length(S4Vectors::subjectHits(hitsEnd2))>0) ){
  				hitsEndApply=tapply(
  				  c(noRead1[S4Vectors::subjectHits(hitsEnd1)], 
  				    noRead2[S4Vectors::subjectHits(hitsEnd2)]),
  					c(S4Vectors::queryHits(hitsEnd1), S4Vectors::queryHits(hitsEnd2)),
  					function(tmp) return(length(unique(tmp))) )
  				hitsBegApply=tapply(
  				  c(noRead1[S4Vectors::subjectHits(hitsBeg1)], 
  				    noRead2[S4Vectors::subjectHits(hitsBeg2)]),
  					c(S4Vectors::queryHits(hitsBeg1), S4Vectors::queryHits(hitsBeg2)), 
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
  		strRead1=sam1[,"strand"]
  		chrRead2=sam2[,"rName"]
  		locRead2=sam2[,"pos"]
  		strRead2=sam2[,"strand"]
  		
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
  				frqEventTmpM, begEventTmpM, endEventTmpM,
  				strand=rep(strRead1[tmp],length(frqEventTmpM))))
  	
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
  				frqEventTmpM, begEventTmpM, endEventTmpM,
  				strand=rep(strRead2[tmp],length(frqEventTmpM))))
  		})
  
  		noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[1]]))))
  		chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			as.character(unlist(tmp[[2]])))))
  		begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[4]]))))
  		endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  			unlist(tmp[[5]]))))
  		strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[1]]))))
  		chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			as.character(unlist(tmp[[2]])))))
  		begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[4]]))))
  		endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  			unlist(tmp[[5]]))))
  		strRead2=as.vector(unlist(sapply(frqEventRead2, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		noVec=c(noRead1,noRead2)
  		chrVec=c(as.character(chrRead1),
  			as.character(chrRead2))
  		begVec=c(begRead1,begRead2)
  		endVec=c(endRead1,endRead2)
  		
  		if(strandSpecific=="stranded"){
  		  strAdj1<- as.character(strRead1)
  		  strAdj2<- as.character(revStr[strRead2])
  		} else if (strandSpecific=="reverse"){
  		  strAdj2<- as.character(strRead2)
  		  strAdj1<- as.character(revStr[strRead1])
  		}
  
  		ref=reference[which(referenceIntronExon=="intron"),]
  		read1GRanges=GenomicRanges::GRanges( seqnames=chrRead1, 
  		  IRanges::IRanges(start=begRead1, end=endRead1), strand=strAdj1)
  		read2GRanges=GenomicRanges::GRanges( seqnames=chrRead2, 
  		  IRanges::IRanges(start=begRead2, end=endRead2), strand=strAdj2)
  
  		refIntBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			IRanges::IRanges(start=ref[,"begin"]-1, end=ref[,"begin"]-1), 
  			strand=ref[,"strand"])
  		refIntEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			IRanges::IRanges(start=ref[,"end"]+1, end=ref[,"end"]+1),
  			strand=ref[,"strand"])
  
  		hitsBegInt1 <- GenomicRanges::findOverlaps(read1GRanges, 
  			refIntBegGRanges, type= "end", ignore.strand=FALSE)
  		hitsEndInt1 <- GenomicRanges::findOverlaps(read1GRanges, 
  			refIntEndGRanges, type= "start", ignore.strand=FALSE)
  		
  		hitsBegInt2 <- GenomicRanges::findOverlaps(read2GRanges, 
  		  refIntBegGRanges, type= "end", ignore.strand=FALSE)
  		hitsEndInt2 <- GenomicRanges::findOverlaps(read2GRanges, 
  		  refIntEndGRanges, type= "start", ignore.strand=FALSE)
  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter1= GenomicRanges::findOverlaps(read1GRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
  				hitsFilter2= GenomicRanges::findOverlaps(read2GRanges, 
  				  repeatGRanges, type= "any", ignore.strand=TRUE)
  				
  				noFilQuHits<- 
  				  intersect(noRead1[queryHits(hitsFilter1)], 
  				            noRead2[queryHits(hitsFilter2)])
  				
  				filtEndInd1= which(! noRead1[S4Vectors::subjectHits(hitsEndInt1)] %in% 
  				  noFilQuHits)
  				filtBegInd1= which(! noRead1[S4Vectors::subjectHits(hitsBegInt1)] %in% 
  				  noFilQuHits)
  				
  				filtEndInd2= which(! noRead2[S4Vectors::subjectHits(hitsEndInt2)] %in% 
  				  noFilQuHits)
  				filtBegInd2= which(! noRead2[S4Vectors::subjectHits(hitsBegInt2)] %in% 
  				  noFilQuHits)
  				
  				
  				hitsBegInt1=hitsBegInt1[filtBegInd1,]
  				hitsEndInt1=hitsEndInt1[filtEndInd1,]
  				hitsBegInt2=hitsBegInt2[filtBegInd2,]
  				hitsEndInt2=hitsEndInt2[filtEndInd2,]
  			}
  		
  			# Filter reads mapped completely to Either introns or exons
  
  			refRes= rep(0,nrow(ref))
  			if(length(c(S4Vectors::queryHits(hitsBegInt1), 
  				S4Vectors::queryHits(hitsEndInt1)))>0 |
  				length(c(S4Vectors::queryHits(hitsBegInt2), 
  				S4Vectors::queryHits(hitsEndInt2)))>0){
  				hitsApply=tapply(
  					c(noRead1[c(S4Vectors::queryHits(hitsBegInt1), 
  						  S4Vectors::queryHits(hitsEndInt1))], 
  						noRead2[c(S4Vectors::queryHits(hitsBegInt2), 
  						  S4Vectors::queryHits(hitsEndInt2))]), 
  					c(S4Vectors::subjectHits(hitsBegInt1), 
  						S4Vectors::subjectHits(hitsEndInt1),
  						S4Vectors::subjectHits(hitsBegInt2), 
  						S4Vectors::subjectHits(hitsEndInt2)),
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
  	  
  	  locRead1=sam1[,"pos"]
  	  chrRead1=sam1[,"rName"]
  	  strRead1=sam1[,"strand"]
  	  chrRead2=sam2[,"rName"]
  	  locRead2=sam2[,"pos"]
  	  strRead2=sam2[,"strand"]
  	  
  	  
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
  	              readNoM=rep(tmp,length(frqEventTmpM)),
  	              strand=rep(strRead1[tmp],length(frqEventTmpN)) ))
  	    
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
  	              readNoM=rep(tmp,length(frqEventTmpM)),
  	              strand=rep(strRead2[tmp],length(frqEventTmpN) )))
  	    
  	  })
  	  
  	  noRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[1]]))))
  	  chrRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[2]]))))
  	  begRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[4]]))))
  	  endRead1= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    unlist(tmp[[5]]))))
  	  strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    as.character(unlist(tmp[[10]])))))
  
  	  noRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[1]]))))
  	  chrRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[2]]))))
  	  begRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[4]]))))
  	  endRead2= as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    unlist(tmp[[5]]))))
  	  strRead2=as.vector(unlist(sapply(frqEventRead2, function(tmp)
  	    as.character(unlist(tmp[[10]])))))
  	  
  	  if(strandSpecific=="stranded"){
  	    strAdj1<- as.character(strRead1)
  	    strAdj2<- as.character(revStr[strRead2])
  	  } else if (strandSpecific=="reverse"){
  	    strAdj2<- as.character(strRead2)
  	    strAdj1<- as.character(revStr[strRead1])
  	  }
  	  
  	  # analyzing reads
  	  noVec=c(noRead1,noRead2)
  	  chrVec=c(as.character(chrRead1),
  	           as.character(chrRead2))
  	  begVec=c(begRead1,begRead2)
  	  endVec=c(endRead1,endRead2)
  	  strVec<- c(strAdj1,strAdj2)
  	  
  	  ref=reference[which(referenceIntronExon=="exon"),]
  	  read1GRanges=GenomicRanges::GRanges( seqnames=chrRead1, 
  	    IRanges::IRanges(start=begRead1, end=endRead1), strand=strAdj1)
  	  read2GRanges=GenomicRanges::GRanges( seqnames=chrRead2, 
  	    IRanges::IRanges(start=begRead2, end=endRead2), strand=strAdj2)
  	  
  	  refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  	    IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"]))
  	  
  	  hitsAllEx1 <- GenomicRanges::findOverlaps(refGRanges, 
  	                                           read1GRanges, type= "within",
  	                                           ignore.strand= FALSE)
  	  hitsAllEx2 <- GenomicRanges::findOverlaps(refGRanges, 
  	                                            read2GRanges, type= "within",
  	                                            ignore.strand= FALSE)
  	  

  	  
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
  	    
  	    readGRangesM1=GenomicRanges::GRanges( seqnames=chrRead1M, 
  	      IRanges::IRanges(start=begRead1M, end=endRead1M))
  	    repeatGRanges= GenomicRanges::GRanges( 
  	      seqnames=repeatsTableToFilter[,"chr"], 
  	      IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  	                       end=repeatsTableToFilter[,"end"], 
  	                       width=repeatsTableToFilter[,"end"]-
  	                         repeatsTableToFilter[,"begin"]+1))
  	    hitsFilter1= GenomicRanges::findOverlaps(readGRangesM1, 
  	      repeatGRanges, type= "any", ignore.strand=TRUE)
  	    hitsFilter2= GenomicRanges::findOverlaps(readGRangesM2, 
  	      repeatGRanges, type= "any", ignore.strand=TRUE)

  	    # noVecFil<- unique(noVecM[queryHits(hitsFilter)])
  	    # filthitsAllExInd<- which(!subjectHits(hitsAllEx) %in% 
  	    #                            subjectHits(hitsAllEx)[
  	    #                              which(noVec[subjectHits(hitsAllEx)]
  	    #                                    %in% noVecFil)])
  	    # hitsAllEx=hitsAllEx[filthitsAllExInd,]
  	    noFilQuHits<- 
  	      intersect(noRead1[queryHits(hitsFilter1)], 
  	                noRead2[queryHits(hitsFilter2)])
     
  	    filtAllInd1= which(! noRead1[S4Vectors::subjectHits(hitsAllEx1)] %in% 
  	                         noFilQuHits)
  	    filtAllInd2= which(! noRead1[S4Vectors::subjectHits(hitsAllEx2)] %in% 
  	                         noFilQuHits)
  	    

  	    
  	    hitsAllEx1=hitsAllEx1[filtAllInd1,]
  	    hitsAllEx2=hitsAllEx2[filtAllInd2,]
  	    # hitsBegInt2=hitsBegInt2[filtBegInd2,]
  	    # hitsEndInt2=hitsEndInt2[filtEndInd2,]
  	    
  	  }
  	  
  	  # Filter reads mapped completely to Either introns or exons
  	  
  	  refRes= rep(0,nrow(ref))
  	  if(length(c(S4Vectors::subjectHits(hitsAllEx1),
  	              S4Vectors::subjectHits(hitsAllEx2)))>0){
  	    hitsApply=tapply(
  	      c(noRead1[S4Vectors::subjectHits(hitsAllEx1)],
  	        noRead2[S4Vectors::subjectHits(hitsAllEx2)]), 
  	      c(S4Vectors::queryHits(hitsAllEx1),
  	        S4Vectors::queryHits(hitsAllEx2)),
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
	
	if(length(intRetRes)==0)
	  intRetRes<- rep(0,nrow(reference))

	if(length(exExRes)==0)
	  exExRes<- rep(0,nrow(reference))
	  
	if(length(intSpanRes)==0)
	  intSpanRes<- rep(0,nrow(reference))

	if(length(exSkipRes)==0)
	  exSkipRes<- rep(0,nrow(reference))

#	if(length(intRetRes)>0)
	if("IntRet" %in% method)
		finalRes<-c(finalRes, intRetRes)
	rm("intRetRes")
	
#	if(length(exExRes)>0){
	if("ExEx" %in% method)
		finalRes<-c(finalRes, exExRes)
	rm("exExRes")

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
interestIntExAnalyseSingleStranded <- function(
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
	strandSpecific,
	excludeFusionReads)
{
  revStr<- c("-", "+", "*")
  names(revStr)<- c("+", "-", "*")
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
	
	sam2= NA
	sam1Flag<- mcols(unlist(readTmpFil))$flag
	sam1FlagMat<-bamFlagAsBitMatrix(sam1Flag)
	readMate<- rep(1,nrow(sam1FlagMat))
	readMate[which(sam1FlagMat[,"isSecondMateRead"]==1)]<- 2

	if(length(limitRanges)>0 & length(readTmpFil)>0){
	  
	  readTmpFirst<- readTmpFil 
	  if(strandSpecific=="stranded"){
	    if(length(which(readMate==1))>0)
	      GenomicRanges::strand(readTmpFirst)[which(readMate==1)]<- 
	        as.character(GenomicRanges::strand(readTmpFirst)[which(readMate==1)])
	    if(length(which(readMate==2))>0)
	      GenomicRanges::strand(readTmpFirst)[which(readMate==2)]<- 
	        as.character(revStr[as.character(GenomicRanges::strand(readTmpFirst))][which(readMate==2)])
	  } else if (strandSpecific=="reverse"){
	    if(length(which(readMate==2))>0)
	      GenomicRanges::strand(readTmpFirst)[which(readMate==2)]<- 
	        as.character(GenomicRanges::strand(readTmpFirst)[which(readMate==2)])
	    if(length(which(readMate==1))>0)
	      GenomicRanges::strand(readTmpFirst)[which(readMate==1)]<- 
	        as.character(revStr[as.character(GenomicRanges::strand(readTmpFirst))][which(readMate==1)])
	  }
	  
	  r1Map<- GenomicRanges::findOverlaps(
	    readTmpFirst, 
	    limitRanges,
	    type= "within",
	    ignore.strand= FALSE)
	  
	  filLimInd<- sort(unique(S4Vectors::queryHits(r1Map)), decreasing=FALSE)
	  readTmpFil<- readTmpFil[filLimInd,]
	}


	
#Extract information of the mapped pieces of the reads
	sam1<- as.data.frame(readTmpFil)[, c("qname", "rname", "strand", 
		"pos", "qual", "cigar")]
	colnames(sam1)<- c("qName","rName","strand","pos","qual",
		"cigar")

	### Incorporate FLAGS to figure what reads are fisrt or second from mate in singletons
	### then adjust strands accordingly.
	### Use flags<-bamFlagAsBitMatrix(mcols(unlist(readTmpFil))$flag)
	### Column 7 and 8 is first or second mate
	### Update what in scanbamwhat with scParam3=Rsamtools::ScanBamParam(
	###what=Rsamtools::scanBamWhat()[c(1,
	###                                3,5,8,13,9, 10, 6, 4, 14, 15, 2)],

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
  		strRead1=sam1[,"strand"]
  		strRead2=0

  		if(strandSpecific=="stranded"){
  		  if(length(which(readMate==2))>0)
  		    strRead1[which(readMate==2)]<- 
  		      as.character(revStr[which(readMate==2)]) 
  		} else if (strandSpecific=="reverse"){
  		  strRead1[which(readMate==1)]<- 
  		    as.character(revStr[strRead1[which(readMate==1)]])
  		}  

  
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
  				frqEventTmpM, begEventTmpM, endEventTmpM,
  				strand=rep(strRead1[tmp],length(frqEventTmpM)),
  				mate=rep(readMate[tmp],length(frqEventTmpM))))
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
  		strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		mateRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[7]])))))
  		noRead2= c()
  		chrRead2= c()
  		begRead2= c()
  		endRead2= c()
  		strRead2= c()
  		strRead2= c()
  		strAdj1<- strRead1
  		
  		##### HERE 1 !!!!
  		strAdj2<- c()
  # analyzing reads
  		noVec=c(noRead1,noRead2)
  		chrVec=c(as.character(chrRead1),
  			as.character(chrRead2))
  		begVec=c(begRead1,begRead2)
  		endVec=c(endRead1,endRead2)
  		strVec=c(as.character(strAdj1),
  		         as.character(strAdj2))
  		exExRes<- c()
  		intRetRes<- c()
  		methodNo= which(method=="ExEx")
  		if(length(methodNo)>0){
  			ref= reference[which(referenceIntronExon=="exon"),]
  			readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  				IRanges::IRanges(start=begVec, end=endVec), strand=strVec)
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"], 
  					width=ref[,"end"]-ref[,"begin"]+1), strand=ref[,"strand"])
  			hits1 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "start", ignore.strand=FALSE)
  			if(logFile!="")
  			  cat( "Part8 ",qNam, "\n",
  			       file=logFile, append=appendLogFile)
  			hits2 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "end", ignore.strand=FALSE)
  			if(logFile!="")
  			  cat( "Part9 ",qNam, "\n", 
  			       file=logFile, append=appendLogFile)
  			hits3 <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "equal", ignore.strand=FALSE)
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
  					repeatGRanges, type= "any", ignore.strand=TRUE)
  	
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
  					refGRanges, type= "within", ignore.strand=FALSE)
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
  				IRanges::IRanges(start=begVec, end=endVec), strand=strVec)
  			refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges(start=ref[,"begin"], end=ref[,"end"], 
  					width=ref[,"end"]-ref[,"begin"]+1), strand=ref[,"strand"])
  			hits <- GenomicRanges::findOverlaps(readGRanges, refGRanges, 
  				type= "any", ignore.strand=FALSE)
  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
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
  				IRanges::IRanges(start=begVec, end=endVec), strand=strVec)
  			refEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"end"], end=ref[,"end"]),
  				strand=ref[,"strand"])
  			refBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  				IRanges::IRanges(start=ref[,"begin"], end=ref[,"begin"]),
  				strand=ref[,"strand"])
  			hitsEnd <- GenomicRanges::findOverlaps(refEndGRanges, readGRanges,
  				type= "within", ignore.strand=FALSE)
  			hitsBeg <- GenomicRanges::findOverlaps(refBegGRanges, readGRanges,
  				type= "within", ignore.strand=FALSE)
  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
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
  	}
    
### HERE !!!
    
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
  		lenRead2=0
  		strRead1=sam1[,"strand"]
  		strRead2=0
  		
  		if(strandSpecific=="stranded"){
  		  if(length(which(readMate==2))>0)
  		    strRead1[which(readMate==2)]<- 
  		      as.character(revStr[which(readMate==2)]) 
  		} else if (strandSpecific=="reverse"){
  		  strRead1[which(readMate==1)]<- 
  		    as.character(revStr[strRead1[which(readMate==1)]])
  		}  

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
  				frqEventTmpM, begEventTmpM, endEventTmpM, 
  				strand=rep(strRead1[tmp],length(frqEventTmpM)),
  				mate=rep(readMate[tmp],length(frqEventTmpM)))) 
  
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
  		strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[6]])))))
  		mateRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  		  as.character(unlist(tmp[[7]])))))
  		noRead2= c()
  		chrRead2= c()
  		begRead2= c()
  		endRead2= c()
  		strRead2= c()
  		strAdj1<- strRead1
  		strAdj2<- c()

  # analyzing reads
  		noVec=c(noRead1,noRead2)
  		chrVec=c(as.character(chrRead1),
  			as.character(chrRead2))
  		begVec=c(begRead1,begRead2)
  		endVec=c(endRead1,endRead2)
  		strVec=c(as.character(strAdj1),
  		         as.character(strAdj2))
  		ref=reference[which(referenceIntronExon=="intron"),]
  		readGRanges=GenomicRanges::GRanges( seqnames=chrVec, 
  			IRanges::IRanges(start=begVec, end=endVec), strand=strVec)
  
  		refIntBegGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			IRanges::IRanges(start=ref[,"begin"]-1, end=ref[,"begin"]-1), 
  			strand=ref[,"strand"])
  		refIntEndGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"], 
  			IRanges::IRanges(start=ref[,"end"]+1, end=ref[,"end"]+1), 
  			strand=ref[,"strand"])
  
  		hitsBegInt <- GenomicRanges::findOverlaps(readGRanges, 
  			refIntBegGRanges, type= "end", ignore.strand=FALSE)
  		hitsEndInt <- GenomicRanges::findOverlaps(readGRanges, 
  			refIntEndGRanges, type= "start", ignore.strand=FALSE)
  			# Filter reads mapped to repeats regions if requested
  			if(length(repeatsTableToFilter)!=0){
  				repeatGRanges= GenomicRanges::GRanges( 
  					seqnames=repeatsTableToFilter[,"chr"], 
  					IRanges::IRanges(start=repeatsTableToFilter[,"begin"],
  						end=repeatsTableToFilter[,"end"], 
  						width=repeatsTableToFilter[,"end"]-
  							repeatsTableToFilter[,"begin"]+1))
  				hitsFilter= GenomicRanges::findOverlaps(readGRanges, 
  					repeatGRanges, type= "any", ignore.strand=TRUE)
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
  	  strRead1=sam1[,"strand"]
  	  chrRead2= c()
  	  strRead2=c()
  	  locRead2=c()
  	  
  	  if(strandSpecific=="stranded"){
  	    if(length(which(readMate==2))>0)
  	      strRead1[which(readMate==2)]<- 
  	        as.character(revStr[which(readMate==2)]) 
  	  } else if (strandSpecific=="reverse"){
  	    strRead1[which(readMate==1)]<- 
  	      as.character(revStr[strRead1[which(readMate==1)]])
  	  }  
  	  
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
  	                 readNoM=rep(tmp,length(frqEventTmpM)),
  	                 strand=rep(strRead1[tmp],length(frqEventTmpN)),
  	                 mate=rep(readMate[tmp],length(frqEventTmpN)),
  	                 chrM=rep(chrRead1[tmp],length(frqEventTmpM))))
  	    
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
  	  strRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    as.character(unlist(tmp[[10]])))))
  	  mateRead1=as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	    as.character(unlist(tmp[[11]])))))
  	  noRead2= c()
  	  chrRead2= c()
  	  begRead2= c()
  	  endRead2= c()
  	  strRead2=c()
  	  strAdj1<- strRead1
  	  
  	  # analyzing reads
  	  noVec=c(noRead1,noRead2)
  	  chrVec=c(as.character(chrRead1),
  	           as.character(chrRead2))
  	  begVec=c(begRead1,begRead2)
  	  endVec=c(endRead1,endRead2)
  	  strVec=c(as.character(strAdj1),
  	           as.character(strAdj2))
  	  
  	  ref=reference[which(referenceIntronExon=="exon"),]
  	  readGRanges=GenomicRanges::GRanges( seqnames=chrVec,
  	    IRanges::IRanges(start=begVec, end=endVec), strand=strVec)
  	  
  	  refGRanges= GenomicRanges::GRanges( seqnames=ref[,"chr"],
        IRanges::IRanges(start=ref[,"begin"], end=ref[,"end"]),
        strand=ref[,"strand"])
  
  	  hitsAllEx <- GenomicRanges::findOverlaps(refGRanges, 
        readGRanges, type= "within", ignore.strand=FALSE)
  
  	  # Filter reads mapped to repeats regions if requested
  	  if(length(repeatsTableToFilter)!=0){
  	    noRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[9]]))))
  	    chrRead1M= as.vector(unlist(sapply(frqEventRead1, function(tmp)
  	      unlist(tmp[[12]]))))
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
          repeatGRanges, type= "any", ignore.strand=TRUE)
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
	
	if(length(intRetRes)==0)
	  intRetRes<- rep(0,nrow(reference))
	
	if(length(exExRes)==0)
	  exExRes<- rep(0,nrow(reference))
	
	if(length(intSpanRes)==0)
	  intSpanRes<- rep(0,nrow(reference))
	
	if(length(exSkipRes)==0)
	  exSkipRes<- rep(0,nrow(reference))
	
	#	if(length(intRetRes)>0)
	if("IntRet" %in% method)
	  finalRes<-c(finalRes, intRetRes)
	rm("intRetRes")
	
	#	if(length(exExRes)>0){
	if("ExEx" %in% method)
	  finalRes<-c(finalRes, exExRes)
	rm("exExRes")
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
	

