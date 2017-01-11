interestAnalyse <-
function(
	reference=u12,
	outDir,
	logFile="",
	pairFiles,
	singleFiles,
	inLoc,
	method=c("IntRet","ExEx"),
	referenceIntronExon=u12[,"int_ex"],
	clusterNo="",
	appendLogFile=TRUE,
	repeatsTableToFilter=c(),
	junctionReadsOnly=FALSE)
{

#Function running on all parallel cores 

	interestIntExAnalyse <- function(){
		junk <- 0 
		done <- 0 
		while (done != 1) {
# Signal being ready to receive a new task 
		Rmpi::mpi.send.Robj(junk,0,1) 

# Receive a task 
		task <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag())
		task_info <- Rmpi::mpi.get.sourcetag() 
		tag <- task_info[2] 
#Begin of parallel process
		if (tag == 1) {
		#library(GenomicRanges)
			fileNumber <- task$fileNumber

	#Define functions
			subjectHits<-  S4Vectors::subjectHits
			queryHits<- S4Vectors::queryHits
			GRanges<- GenomicRanges::GRanges
			findOverlaps<- GenomicRanges::findOverlaps
			IRanges<- IRanges::IRanges
	
			if(logFile!="")
				cat( "InERESt:interestAnalyse: Cluster number ",
					Rmpi::mpi.comm.rank(), " analyzing file ",
					inFiles[fileNumber], "\n",file=logFile, append=TRUE)
			cat( "InERESt:interestAnalyse: Cluster number ",
				Rmpi::mpi.comm.rank(), " analyzing file ", 
				inFiles[fileNumber], "\n")

			count=Rmpi::mpi.comm.rank()
			if(paired[fileNumber]){
				sam1=read.table(
					paste(inLoc,"/read1/",inFiles[fileNumber], sep=""),
					header=FALSE, stringsAsFactors=FALSE, comment.char="")
				colnames(sam1)= c("qName","rName","strand","pos","qual",
					"cigar", "rNext","posNext")
				sam2=read.table(
					paste(inLoc,"/read2/",inFiles[fileNumber], sep=""),
					header=FALSE, stringsAsFactors=FALSE, comment.char="")
				colnames(sam2)=c("qName", "rName", "strand", "pos", "qual",
					"cigar", "rNext","posNext")

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
				chrRead2=sam1[,"rNext"]
				locRead2=sam2[,"pos"]

				frqEventRead1=lapply(1:length(sam1$cigar), function (tmp) {
					frqEventTmp= as.numeric(unlist(strsplit(sam1$cigar[tmp],
						split="[A-Z]")))
					names(frqEventTmp)= unlist(strsplit(sam1$cigar[tmp], 
						split=""))[grep('[A-Z]', unlist(strsplit(sam1$cigar[tmp]
							,split="")))]
					frqEventTmp=frqEventTmp[!is.na(match(names(frqEventTmp),
						c("N","D","M")))]
					begEventTmp=rep(NA,length(frqEventTmp))
					endEventTmp=rep(NA,length(frqEventTmp))
					begEventTmp[1]=locRead1[tmp]
					endEventTmp[1]=begEventTmp[1]+frqEventTmp[1]-1

					if(length(frqEventTmp)>1) 
						for(cnt in 2:length(frqEventTmp)){
							begEventTmp[cnt]=endEventTmp[cnt-1]+1
							endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
						}
					frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					list(readNo=rep(tmp,length(frqEventTmpM)), 
						chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
						frqEventTmpM, begEventTmpM, endEventTmpM)
	
				})
	
				frqEventRead2=lapply(1:length(sam2$cigar), function (tmp) {
					frqEventTmp= as.numeric(unlist(strsplit(sam2$cigar[tmp],
						split="[A-Z]")))
					names(frqEventTmp)= unlist(strsplit(sam2$cigar[tmp], 
						split=""))[grep('[A-Z]', unlist(strsplit(sam2$cigar[tmp]
						, split="")))]
					frqEventTmp= frqEventTmp[!is.na(match(names(frqEventTmp),
						c("N","D","M")))]
					begEventTmp=rep(NA,length(frqEventTmp))
					endEventTmp=rep(NA,length(frqEventTmp))
					begEventTmp[1]=locRead2[tmp]
					endEventTmp[1]=begEventTmp[1]+frqEventTmp[1]-1
					if(length(frqEventTmp)>1) 
						for(cnt in 2:length(frqEventTmp)){
							begEventTmp[cnt]=endEventTmp[cnt-1]+1
							endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
						}
					frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp),
						"M"))]
					list(readNo=rep(tmp,length(frqEventTmpM)), 
						chr=rep(chrRead2[tmp],length(frqEventTmpM)), 
						frqEventTmpM, begEventTmpM, endEventTmpM)
				})

				noRead1= unlist(sapply(frqEventRead1, function(tmp)
					unlist(tmp[[1]])))
				chrRead1= unlist(sapply(frqEventRead1, function(tmp)
					unlist(tmp[[2]])))
				begRead1= unlist(sapply(frqEventRead1, function(tmp)
					unlist(tmp[[4]])))
				endRead1= unlist(sapply(frqEventRead1, function(tmp)
					unlist(tmp[[5]])))
				noRead2= unlist(sapply(frqEventRead2, function(tmp)
					unlist(tmp[[1]])))
				chrRead2= unlist(sapply(frqEventRead2, function(tmp)
					unlist(tmp[[2]])))
				begRead2= unlist(sapply(frqEventRead2, function(tmp)
					unlist(tmp[[4]])))
				endRead2= unlist(sapply(frqEventRead2, function(tmp)
					unlist(tmp[[5]])))
	
			} else {
	
				sam1= read.table(
					paste(inLoc,"/single/",inFiles[fileNumber],sep=""), 
					header=FALSE, stringsAsFactors=FALSE, comment.char="")
		
				colnames(sam1)= c("qName","rName","strand","pos","qual","cigar"
					,"rNext","posNext")
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
					begEventTmp=rep(NA,length(frqEventTmp))
					endEventTmp=rep(NA,length(frqEventTmp))
					begEventTmp[1]= locRead1[tmp]
					endEventTmp[1]= begEventTmp[1]+frqEventTmp[1]-1
					if(length(frqEventTmp)>1) 
						for(cnt in 2:length(frqEventTmp)){
							begEventTmp[cnt]=endEventTmp[cnt-1]+1
							endEventTmp[cnt]=begEventTmp[cnt]+frqEventTmp[cnt]-1	
						}
					frqEventTmpM=frqEventTmp[!is.na(match(names(frqEventTmp), 
						"M"))]
					begEventTmpM=begEventTmp[!is.na(match(names(frqEventTmp), 
						"M"))]
					endEventTmpM=endEventTmp[!is.na(match(names(frqEventTmp), 
						"M"))]
					list(readNo=rep(tmp,length(frqEventTmpM)), 
						chr=rep(chrRead1[tmp],length(frqEventTmpM)), 
						frqEventTmpM, begEventTmpM, endEventTmpM)
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
			
			}
			noVec=c(noRead1,noRead2)
			chrVec=c(chrRead1,chrRead2)
			begVec=c(begRead1,begRead2)
			endVec=c(endRead1,endRead2)
			methodNo=which(method=="ExEx")
			if(length(methodNo)>0){
				tmpLoc=paste(outDir,method[methodNo],sep="/")
				dir.create(tmpLoc)
				ref=reference[referenceIntronExon=="exon",]
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
				write.table(matrix(refRes, ncol=1) , 
					paste(tmpLoc,"/","tmpSam","_",fileNumber,'.frq', sep=''), 
					append=FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE,
					sep='\t' )
			}
			methodNo=which(method=="IntRet")
			if(length(methodNo)>0){
				tmpLoc=paste(outDir,method[methodNo],sep="/")
				if(!dir.exists(tmpLoc))
					dir.create(tmpLoc)
				ref=reference[referenceIntronExon=="intron",]
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
				write.table(matrix(refRes, ncol=1) , 
					paste(tmpLoc,"/","tmpSam","_",fileNumber,'.frq', sep=''), 
					append=FALSE, col.names=FALSE, 
					row.names=FALSE, quote=FALSE, sep='\t' )
			}

		} # End of if tag==1
		else if (tag == 2) {
			done <- 1
		}
	
		} #End of while
	
		Rmpi::mpi.send.Robj(junk,0,3)
	    
	}

# End of function running on all parallel cores

# Finalizing function
	.Last <- function(){
	    if (is.loaded("mpi_initialize")){
        	if (Rmpi::mpi.comm.size(1) > 0){
				Rmpi::mpi.close.Rslaves(1)
        	}
	    }
	}
# End of finalizing  function

	if(logFile!="")
		cat( "InERESt:interestAnalyse: Begins ...\n", file=logFile, 
			append=appendLogFile)
	cat( "InERESt:interestAnalyse: Begins ...\n")

	if(as.character(class(reference))=="GRanges"){
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
				end=GenomicRanges::end(reference), 
				strand=as.character(GenomicRanges::strand(reference)))
		}
		reference=tmpReference
	}


	time1=Sys.time()
	inFiles=c(pairFiles,singleFiles)
	paired=c(rep(TRUE,length(pairFiles)), rep(FALSE,length(singleFiles)))

	# Opening the parallel computer cores
	if(Rmpi::mpi.comm.size()>0){
		Rmpi::mpi.close.Rslaves()
	}
	if(clusterNo!=""){ 
		Rmpi::mpi.spawn.Rslaves(nslaves=clusterNo)
	} else {
		clusterNo=Rmpi::mpi.universe.size()
		Rmpi::mpi.spawn.Rslaves(nslaves=clusterNo)
	}


	if (Rmpi::mpi.comm.size() < 2) {
		msg= 
"InERESt:interestAnalyse: Error: More slave processes are required.\n"
		if(logFile!="")
			cat( msg, file=logFile, append=TRUE)
		cat( msg, file=logFile, append=TRUE)
		Rmpi::mpi.quit()
	}

# Seding required objects to the slaves

	Rmpi::mpi.bcast.Robj2slave(reference)
	Rmpi::mpi.bcast.Robj2slave(logFile)
	Rmpi::mpi.bcast.Robj2slave(inFiles)
	Rmpi::mpi.bcast.Robj2slave(outDir)
	Rmpi::mpi.bcast.Robj2slave(inLoc)
	Rmpi::mpi.bcast.Robj2slave(paired)
	Rmpi::mpi.bcast.Robj2slave(method)
	Rmpi::mpi.bcast.Robj2slave(referenceIntronExon)
	Rmpi::mpi.bcast.Robj2slave(repeatsTableToFilter)
	Rmpi::mpi.bcast.Robj2slave(junctionReadsOnly)

	if(length(which(method=="ExEx" | method=="IntRet" ))>0){
		junk <- 0 
		closed_slaves <- 0 
		n_slaves <- Rmpi::mpi.comm.size()-1
		if(logFile!="")
			cat( "InERESt:interestAnalyse: Read counting begins ...\n", 
				file=logFile, append=TRUE)
		cat( "InERESt:interestAnalyse: Read counting begins ...\n")
		#Sending required functions
		Rmpi::mpi.bcast.Robj2slave(interestIntExAnalyse)
		#Analyze slaves
		Rmpi::mpi.bcast.cmd(interestIntExAnalyse())

		# Create task list
		tasks <- vector('list')
		for (i in 1:length(inFiles)) {
		    tasks[[i]] <- list(fileNumber=i)
		}


		#processing the tasks

		while (closed_slaves < n_slaves) { 
			# Receive a message from a slave 
			message <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),
				Rmpi::mpi.any.tag()) 
			message_info <- Rmpi::mpi.get.sourcetag() 
			slave_id <- message_info[1] 
			tag <- message_info[2] 

			if (tag == 1) { 
				# slave is ready for a task. Give it the next task, or tell it  
				# tasks are done if there are none. 
        			if (length(tasks) > 0) { 
					# Send a task, and then remove it from the task list 
					Rmpi::mpi.send.Robj(tasks[[1]], slave_id, 1); 
					tasks[[1]] <- NULL 
				}
				else { 
					Rmpi::mpi.send.Robj(junk, slave_id, 2) 
				}
			}
			else if (tag == 3) { 
				# A slave has closed down. 
				closed_slaves <- closed_slaves + 1 
			} 
		} 
		time2=Sys.time()
		runTime=difftime(time2,time1, units="secs")
		if(logFile!="")
			cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
				runTime," secs\n", file=logFile, append=TRUE)
		cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
			runTime," secs\n")


	}


	if(logFile!="")
		cat( "InERESt:interestAnalyse: Finalizing!\n", file=logFile, 
			append=TRUE)
	cat( "InERESt:interestAnalyse: Finalizing!\n", file=logFile, append=TRUE)

	.Last()
}
