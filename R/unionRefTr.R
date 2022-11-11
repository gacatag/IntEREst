unionRefTr<- function( referenceChr, 
    referenceBegin, 
    referenceEnd, 
	referenceTr,
	referenceIntronExon,
	intronExon="exon", 
	silent=FALSE){

	loc<- paste(referenceChr, referenceBegin, referenceEnd, sep='/')
	ind<- which(referenceIntronExon %in% intronExon)
	tr<- as.character(referenceTr)
	if(any((intronExon=="exon"|intronExon=="intron") & 
		length(intronExon)==1)){
		tr<- tr[ind]
		loc<-loc[ind]
		referenceIntronExon<- referenceIntronExon[ind]
	}
	locTrMap<- tapply(tr,loc,unique)
	uniLoc<- names(locTrMap)
	tmpMat<- matrix(unlist(strsplit(uniLoc, split="/")), ncol=3, 
		byrow=TRUE)
	res<-cbind( tmpMat, as.vector(unlist(sapply(locTrMap, function(tmp) 
		paste(tmp, collapse=",")))) )
	colnames(res)<- c("chr", "begin", "end", "transcripts_id")

	genList<- strsplit(res[,"transcripts_id"], split=',')

	genVec<- unlist(genList)
	genUni<-unique(genVec)	
	indVecOrig<- unlist(lapply(1:length(genList), function(tmp){ 
		rep(tmp, length(unlist(genList[tmp]))) 
		} ))
	indVec<-indVecOrig
	lenGenUni<-length(genUni)
	lapplyRes<- lapply(1:lenGenUni, function(i){
		if(i %% 1000 == 1 & (!silent)) 
			print(paste(i, lenGenUni, sep="/"))
		genTmp<-genUni[i]
		if(length(genTmp)>0){
			matchGen<-which(genVec==genTmp)
			indTmp<-indVec[matchGen]
			if(length(indTmp)>0 & length(which(is.na(indTmp)))==0){
				matchInds<-which(indVec%in%indTmp)
				indVec[matchInds]<<-min(indTmp)
			}
		}
	})
	genAnnoTmp<- as.character(
		tapply(genVec, indVec, function(tmp) 
			paste(unique(tmp), collapse=",")))
	if(length(which(intronExon%in%"intron"))>0){
		intExRes<- tapply(referenceIntronExon, loc, head, 1)
	} else if(length(intronExon)==1) {
		intExRes<-as.character(rep(intronExon, nrow(res)))
	}
	names(genAnnoTmp)<- unique(indVec)
	names(intExRes)<- unique(indVec)

	resDf<- data.frame(chr= as.character(res[,"chr"]),
		begin= as.numeric(res[,"begin"]),
		end= as.numeric(res[,"end"]),
		transcripts_id= as.character(
			tapply(indVec, indVecOrig, function(tmp) 
				unique(genAnnoTmp[as.character(tmp)]))),
		int_ex=intExRes)
	resDf<- resDf[order(as.character(resDf[,"transcripts_id"]), 
		decreasing=FALSE),]
	return(resDf)
}



