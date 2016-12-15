interestAnalyse.sequential <-
function(
reference=u12,
outDir,
logFile,
pairFiles,
singleFiles,
inLoc,
method=c("IntRet","ExEx"),
referenceIntronExon=u12[,"int_ex"],
appendLogFile=TRUE,
repeatsTableToFilter=c(),
junctionReadsOnly
)
{
	#Defining the required functions
	subjectHits<-  S4Vectors::subjectHits
	queryHits<- S4Vectors::queryHits
	GRanges<- GenomicRanges::GRanges
	findOverlaps<- GenomicRanges::findOverlaps
	IRanges<- IRanges::IRanges

	if(as.character(class(reference))=="GRanges"){
		if(length(reference@ranges@NAMES)){
			tmpReference=data.frame(chr=as.character(reference@seqnames), begin=reference@ranges@start, end=reference@ranges@start+reference@ranges@width-1, strand=reference@strand, 
			names=reference@ranges@NAMES)
		} else{
			tmpReference=data.frame(chr=as.character(reference@seqnames), begin=reference@ranges@start, end=reference@ranges@start+reference@ranges@width-1, strand=reference@strand)
		}
		reference=tmpReference
	}

	cat( "InERESt:interestAnalyse.sequential: Read counting begins ...\n", file=logFile, append=appendLogFile)
	time1=Sys.time()
	inFiles=c(pairFiles,singleFiles)
	paired=c(rep(TRUE,length(pairFiles)), rep(FALSE,length(singleFiles)))


	for (i in 1:length(inFiles)) {
		fileNumber=i
		cat( "InERESt:interestAnalyse.sequential: Analyzing paired/single files ",i," / ", length(inFiles), "\n", file=logFile, append=TRUE)
		if(paired[fileNumber]){
			sam1=read.table(paste(inLoc,"/read1/",inFiles[fileNumber], sep=""), header=FALSE, stringsAsFactors=FALSE, comment.char="")
			colnames(sam1)=c("qName","rName","strand","pos","qual", "cigar", "rNext","posNext")
			sam2=read.table(paste(inLoc,"/read2/",inFiles[fileNumber], sep=""), header=FALSE, stringsAsFactors=FALSE, comment.char="")
			colnames(sam2)=c("qName","rName","strand","pos","qual", "cigar", "rNext","posNext")
			lenRead1= unlist(lapply(sam1$cigar, function(tmp) {
				frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
				names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]', unlist(strsplit(tmp,split="")))]
				return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
			}))
			lenRead2= unlist(lapply(sam2$cigar, function(tmp) {
				frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
				names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]', unlist(strsplit(tmp,split="")))]
				return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))
			}))
			locRead1=sam1[,"pos"]
			chrRead1=sam1[,"rName"]
			chrRead2=sam1[,"rNext"]
			locRead2=sam2[,"pos"]
	
			frqEventRead1=lapply(1:length(sam1$cigar), function (tmp) {frqEventTmp=as.numeric(unlist(strsplit(sam1$cigar[tmp],split="[A-Z]")))
				names(frqEventTmp)=unlist(strsplit(sam1$cigar[tmp], split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp], split="")))]
				frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),c("N","D","M")))]
				begEventTmp=rep(NA,length(frqEventTmp))
				endEventTmp=rep(NA,length(frqEventTmp))
				begEventTmp[1]=locRead1[tmp]
				endEventTmp[1]=begEventTmp[1]+frqEventTmp[1]-1
				if(length(frqEventTmp)>1) 
					for(cnt in 2:length(frqEventTmp)){
						begEventTmp[cnt]=endEventTmp[cnt-1]+1
						endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
					}
				#extract M
				frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				list(readNo=rep(tmp,length(frqEventTmpM)), chr=rep(chrRead1[tmp],length(frqEventTmpM)), frqEventTmpM, begEventTmpM, endEventTmpM)
			})
			frqEventRead2=lapply(1:length(sam2$cigar), function (tmp) {frqEventTmp=as.numeric(unlist(strsplit(sam2$cigar[tmp],split="[A-Z]")))
				names(frqEventTmp)=unlist(strsplit(sam2$cigar[tmp], split=""))[grep('[A-Z]', unlist(strsplit(sam2$cigar[tmp], split="")))]
				frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),c("N","D","M")))]
				begEventTmp=rep(NA,length(frqEventTmp))
				endEventTmp=rep(NA,length(frqEventTmp))
				begEventTmp[1]=locRead2[tmp]
				endEventTmp[1]=begEventTmp[1]+frqEventTmp[1]-1
				if(length(frqEventTmp)>1) 
					for(cnt in 2:length(frqEventTmp)){
						begEventTmp[cnt]=endEventTmp[cnt-1]+1
						endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
					}
				#extract M
				frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				list(readNo=rep(tmp,length(frqEventTmpM)), chr=rep(chrRead2[tmp],length(frqEventTmpM)), frqEventTmpM, begEventTmpM, endEventTmpM)
			})
			noRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[1]])))
			chrRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[2]])))
			begRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[4]])))
			endRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[5]])))
			noRead2=unlist(sapply(frqEventRead2, function(tmp)unlist(tmp[[1]])))
			chrRead2=unlist(sapply(frqEventRead2, function(tmp)unlist(tmp[[2]])))
			begRead2=unlist(sapply(frqEventRead2, function(tmp)unlist(tmp[[4]])))
			endRead2=unlist(sapply(frqEventRead2, function(tmp)unlist(tmp[[5]])))

		} else {

			sam1=read.table(paste(inLoc,"/single/",inFiles[fileNumber],sep=""), header=FALSE, stringsAsFactors=FALSE, comment.char="")
	
			colnames(sam1)=c("qName","rName","strand","pos","qual", "cigar", "rNext","posNext")
			sam2=NA
	
			lenRead1= unlist(lapply(sam1$cigar, function(tmp) {
				frqEvent=as.numeric(unlist(strsplit(tmp,split="[A-Z]")))
				names(frqEvent)=unlist(strsplit(tmp,split=""))[grep('[A-Z]', unlist(strsplit(tmp,split="")))]
				return(sum(frqEvent[grep("[N|M|D]",names(frqEvent))]))		
			}))
			lenRead2=0
			locRead1=sam1[,"pos"]
			chrRead1=sam1[,"rName"]
			chrRead2=0
			locRead2=0
			frqEventRead1=lapply(1:length(sam1$cigar), function (tmp) {frqEventTmp=as.numeric(unlist(strsplit(sam1$cigar[tmp],split="[A-Z]")))
				names(frqEventTmp)=unlist(strsplit(sam1$cigar[tmp], split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp], split="")))]
				frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),c("N","D","M")))]
				begEventTmp=rep(NA,length(frqEventTmp))
				endEventTmp=rep(NA,length(frqEventTmp))
				begEventTmp[1]=locRead1[tmp]
				endEventTmp[1]=begEventTmp[1]+frqEventTmp[1]-1
				if(length(frqEventTmp)>1) 
					for(cnt in 2:length(frqEventTmp)){
						begEventTmp[cnt]=endEventTmp[cnt-1]+1
						endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
					}
				#extract M
				frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp), "M"))]
				list(readNo=rep(tmp,length(frqEventTmpM)), chr=rep(chrRead1[tmp],length(frqEventTmpM)), frqEventTmpM, begEventTmpM, endEventTmpM)
			})

			frqEventRead2=NA
			noRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[1]])))
			chrRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[2]])))
			begRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[4]])))
			endRead1=unlist(sapply(frqEventRead1, function(tmp)unlist(tmp[[5]])))
			noRead2=c()
			chrRead2=c()
			begRead2=c()
			endRead2=c()
		
		}
		noVec=c(noRead1,noRead2)
		chrVec=c(chrRead1,chrRead2)
		begVec=c(begRead1,begRead2)
		endVec=c(endRead1,endRead2)
		methodNo=which(method=="ExEx")
		if(length(methodNo)>0){
			tmpLoc=paste(outDir,method[methodNo],sep="/")
			if(!dir.exists(tmpLoc))
				dir.create(tmpLoc)
			ref=reference[referenceIntronExon=="exon",]
			readGRanges=GRanges( seqnames=chrVec, IRanges(start=begVec, end=endVec))
			refGRanges= GRanges( seqnames=ref[,"chr"], IRanges(start=ref[,"begin"], end=ref[,"end"], width=ref[,"end"]-ref[,"begin"]+1))
			hits1 <- findOverlaps(readGRanges, refGRanges, type= "start")
			hits2 <- findOverlaps(readGRanges, refGRanges, type= "end")
			hits3 <- findOverlaps(readGRanges, refGRanges, type= "equal")
			hitsSubject<- c()
			hitsQuery<- c()
			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GRanges( seqnames=repeatsTableToFilter[,"chr"], IRanges(start=repeatsTableToFilter[,"begin"],
					end=repeatsTableToFilter[,"end"], width=repeatsTableToFilter[,"end"]-repeatsTableToFilter[,"begin"]+1))
				hitsFilter= findOverlaps(readGRanges, repeatGRanges, type= "any")

				filtInd1=which(is.na(match(queryHits(hits1), queryHits(hitsFilter))))
				filtInd2=which(is.na(match(queryHits(hits2), queryHits(hitsFilter))))
				filtInd3=which(is.na(match(queryHits(hits3), queryHits(hitsFilter))))

				hits1=hits1[filtInd1,]
				hits2=hits2[filtInd2,]
				hits3=hits3[filtInd3,]
			}

			# Also include the reads that are fully mapped to exons or introns, if asked by the user
			if(!junctionReadsOnly){
				hitsExtra= findOverlaps(readGRanges, refGRanges, type= "within")

				hitsSubject=subjectHits(hitsExtra)
				hitsQuery=queryHits(hitsExtra)
			}
			hitsSubject=c(hitsSubject, subjectHits(hits1),subjectHits(hits2),subjectHits(hits3))
			hitsQuery=c(hitsQuery, queryHits(hits1),queryHits(hits2),queryHits(hits3))
			refRes= rep(0,nrow(ref))
			hitsApply=tapply(noVec[hitsQuery], hitsSubject, function(tmp) return(length(unique(tmp))) )
			refRes[as.numeric(names(hitsApply))]=as.vector(hitsApply)
			write.table(matrix(refRes, ncol=1) , paste(tmpLoc,"/","tmpSam","_",fileNumber,'.frq', sep=''), append=FALSE, col.names=FALSE, 
				row.names=FALSE, quote=FALSE, sep='\t' )
		}
		methodNo=which(method=="IntRet")
		if(length(methodNo)>0){
			tmpLoc=paste(outDir,method[methodNo],sep="/")
			if(!dir.exists(tmpLoc))
				dir.create(tmpLoc)
			ref=reference[referenceIntronExon=="intron",]
			readGRanges=GRanges( seqnames=chrVec, IRanges(start=begVec, end=endVec))
			refGRanges= GRanges( seqnames=ref[,"chr"], IRanges(start=ref[,"begin"], end=ref[,"end"], 
				width=ref[,"end"]-ref[,"begin"]+1))

			hits <- findOverlaps(readGRanges, refGRanges, type= "any")

			# Filter reads mapped to repeats regions if requested
			if(length(repeatsTableToFilter)!=0){
				repeatGRanges= GRanges( seqnames=repeatsTableToFilter[,"chr"], IRanges(start=repeatsTableToFilter[,"begin"], 
					end=repeatsTableToFilter[,"end"], width=repeatsTableToFilter[,"end"]-repeatsTableToFilter[,"begin"]+1))
				hitsFilter= findOverlaps(readGRanges, repeatGRanges, type= "any")
				filtInd= which(is.na(match(queryHits(hits), queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}

			# Filter reads mapped completely to Either introns or exons
			if(junctionReadsOnly){
				hitsFilter= findOverlaps(readGRanges, refGRanges, type= "within")
				filtInd=which(is.na(match(queryHits(hits), queryHits(hitsFilter))))
				hits=hits[filtInd,]
			}
			refRes= rep(0,nrow(ref))
			hitsApply=tapply(noVec[queryHits(hits)], subjectHits(hits), function(tmp) return(length(unique(tmp))) )
			refRes[as.numeric(names(hitsApply))]=as.vector(hitsApply)
			write.table(matrix(refRes, ncol=1) , paste(tmpLoc,"/","tmpSam","_",fileNumber,'.frq', sep=''), 
				append=FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t' )
		}


	}


	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")

	cat( "InERESt:interestAnalyse.sequential: Read counting ends. Running time: ",runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse.sequential: Finalizing!\n", file=logFile, append=TRUE)
}
